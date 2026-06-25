#!/usr/bin/env python3
"""
Build-time code generator for HPCReact geochemistry reaction systems.
Generates constexpr C++ code from JSON input for GPU efficiency.
"""

import argparse
import json
import os
from typing import List, Dict
from datetime import datetime


class ReactionSystemGenerator:
    def __init__(self, config: Dict):
        self.system_name = config['systemName']
        self.namespace = config['namespace']

        # Primary species (always required)
        self.primary_species = config['primarySpecies']

        # Secondary species - can be provided or auto-inferred
        if 'secondarySpecies' in config:
            # Explicitly provided
            self.secondary_species = config['secondarySpecies']
            self.num_equilibrium_reactions = len(self.secondary_species)
        elif 'numEquilibriumReactions' in config and config['numEquilibriumReactions'] > 0:
            # Auto-infer from first N reactions
            num_eq = config['numEquilibriumReactions']
            print(f"  Info: Auto-inferring {num_eq} secondary species from equilibrium reactions")
            self.secondary_species = []
            self.num_equilibrium_reactions = num_eq
            # Will be populated after parsing reactions
            self._auto_infer_secondary = True
        else:
            # Pure kinetic system (no equilibrium reactions)
            self.secondary_species = []
            self.num_equilibrium_reactions = 0
            self._auto_infer_secondary = False

        # Parse reactions first (before building species list if auto-inferring)
        self.reactions = []
        for rxn_data in config['reactions']:
            if 'equation' in rxn_data:
                # Parse equation format - but we need species list
                # For auto-infer, we'll do this in two passes
                rxn = rxn_data  # Store raw for now
            else:
                # Use stoichiometry array directly
                rxn = rxn_data
            self.reactions.append(rxn)

        # If auto-inferring, extract secondary species from first N reactions
        if hasattr(self, '_auto_infer_secondary') and self._auto_infer_secondary:
            self._infer_secondary_species_from_reactions()

        # Now build combined species list
        if self.secondary_species:
            self.species = self.secondary_species + self.primary_species
        else:
            self.species = self.primary_species
            print("  Info: Pure kinetic system (no secondary species)")

        # Parse equations now that we have the species list
        for i, rxn in enumerate(self.reactions):
            if isinstance(rxn, dict) and 'equation' in rxn and 'stoichiometry' not in rxn:
                self.reactions[i] = self._parse_equation(rxn)

        # Validate: must have one reaction per secondary species + kinetic reactions
        if self.secondary_species and len(self.reactions) < len(self.secondary_species):
            raise ValueError(f"Need at least {len(self.secondary_species)} reactions (one per secondary species), got {len(self.reactions)}")

    def _infer_secondary_species_from_reactions(self):
        """Auto-infer secondary species from equilibrium reactions (first N reactions)"""
        import re

        for i in range(self.num_equilibrium_reactions):
            if i >= len(self.reactions):
                raise ValueError(f"Cannot infer: only {len(self.reactions)} reactions but need {self.num_equilibrium_reactions}")

            rxn = self.reactions[i]
            if 'equation' not in rxn:
                raise ValueError(f"Cannot auto-infer from reaction {i+1}: missing 'equation' field")

            equation = rxn['equation']

            # Extract left side (reactants)
            if '=' not in equation:
                raise ValueError(f"Reaction {i+1}: equation must contain '='")

            left, _ = equation.split('=', 1)
            left = left.strip()

            # Split by + (with spaces) and get terms
            terms = re.split(r'\s+\+\s+', left)

            # Find the secondary species (first non-water, non-solid term on left)
            secondary_species = None
            for term in terms:
                term = term.strip()
                if not term or term in ['H2O', 'H₂O']:
                    continue

                # Remove leading coefficient
                if term[0].isdigit():
                    while term and term[0].isdigit():
                        term = term[1:]
                    term = term.strip()

                # First species on left side that's not water is the secondary species
                # Check if it's in primary species (if so, skip)
                if term not in self.primary_species:
                    secondary_species = term
                    break

            if not secondary_species:
                raise ValueError(f"Reaction {i+1}: Could not identify secondary species on left side of '{equation}'")

            self.secondary_species.append(secondary_species)
            print(f"    Inferred secondary species {i+1}: {secondary_species}")

    def _parse_equation(self, rxn_data: Dict) -> Dict:
        """Parse reaction equation string into stoichiometry"""
        equation = rxn_data['equation']

        # Split into left (reactants) and right (products)
        if '=' not in equation:
            raise ValueError(f"Equation missing '=': {equation}")

        left, right = equation.split('=')
        left = left.strip()
        right = right.strip()

        # Initialize stoichiometry array
        stoich = [0] * len(self.species)

        # Parse left side (negative coefficients)
        self._parse_side(left, stoich, sign=-1)

        # Parse right side (positive coefficients)
        self._parse_side(right, stoich, sign=1)

        return {
            'stoichiometry': stoich,
            'comment': equation,  # Use equation as comment
            'equilibriumConstant': rxn_data['equilibriumConstant'],
            'forwardRate': rxn_data.get('forwardRate', 1.0e10),
            'reverseRate': rxn_data.get('reverseRate', 1.0e10),
            'mobileFlag': rxn_data.get('mobileFlag', 1),
            'isEquilibrium': rxn_data.get('isEquilibrium', True)
        }

    def _parse_side(self, side: str, stoich: List[int], sign: int):
        """Parse one side of equation"""
        import re

        # Skip if empty
        if not side:
            return

        # Split by + surrounded by spaces (to avoid splitting H+ or +2 charges)
        # Match: " + " but not "+" in "H+" or "Ca+2"
        terms = re.split(r'\s+\+\s+', side)

        for term in terms:
            term = term.strip()
            if not term:
                continue

            # Skip water (common spectator species)
            if term in ['H2O', 'H₂O']:
                continue

            # Try to find coefficient only if term starts with digit followed by space or capital letter
            # Examples: "2 Cl-" -> coeff=2, species="Cl-"
            #           "2Cl-" -> coeff=2, species="Cl-"
            #           "Cl-" -> coeff=1, species="Cl-"
            #           "CaCl2" -> coeff=1, species="CaCl2" (no leading digit)

            # First check: does it start with digit(s) followed by space or uppercase?
            if term[0].isdigit():
                # Find where digits end
                i = 0
                while i < len(term) and term[i].isdigit():
                    i += 1
                coeff = int(term[:i])
                species = term[i:].strip()
            else:
                coeff = 1
                species = term

            # Find species index
            if species in self.species:
                idx = self.species.index(species)
                stoich[idx] += sign * coeff
            else:
                # Species not in primary or secondary list
                # Could be a solid/mineral (e.g., CaCO3) or water - skip it
                print(f"  Info: Skipping species '{species}' (likely solid/mineral or solvent)")

    def format_scientific(self, value: float) -> str:
        """Format float in scientific notation: 1.23E+04"""
        if value == 0.0:
            return "0.0"
        s = f"{value:.2E}"
        return s.replace('e', 'E').replace('E+0', 'E+').replace('E-0', 'E-')

    def generate_hpp(self) -> str:
        """Generate complete .hpp file with constexpr arrays"""
        num_reactions = len(self.reactions)
        num_species = len(self.species)

        # Determine number of equilibrium reactions
        if self.num_equilibrium_reactions is not None:
            # Explicitly specified
            num_equilibrium = self.num_equilibrium_reactions
        else:
            # Auto-detect: count reactions with "isEquilibrium": true
            # or fall back to forwardRate >= 1e9
            num_equilibrium = sum(1 for r in self.reactions
                                if r.get('isEquilibrium', False) or
                                   (r.get('forwardRate', 1e10) >= 1e9 and
                                    r.get('reverseRate', 1e10) >= 1e9))

        year = datetime.now().year
        max_species_len = max(len(s) for s in self.species)

        hpp = f"""/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: (BSD-3-Clause)
 *
 * Copyright (c) {year}- Lawrence Livermore National Security LLC
 * All rights reserved
 *
 * See top level LICENSE files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#pragma once

#include "../reactionsSystems/Parameters.hpp"

namespace hpcReact
{{

namespace geochemistry
{{
// turn off uncrustify to allow for better readability of the parameters
// *****UNCRUSTIFY-OFF******

// ################################## {self.system_name} rxn set ##################################
namespace {self.namespace}
{{

constexpr CArrayWrapper<signed char, {num_reactions}, {num_species}> stoichMatrix =
{{ //  """

        # Species header
        for i, sp in enumerate(self.species):
            hpp += f"{sp:>{max_species_len+2}}"
        hpp += "  \n"

        # Stoichiometry matrix rows
        for i, rxn in enumerate(self.reactions):
            hpp += "    {  "
            for j, coeff in enumerate(rxn['stoichiometry']):
                hpp += f"{coeff:>{max_species_len+2}}"
                if j < num_species - 1:
                    hpp += ","
            hpp += "   }"
            if i < num_reactions - 1:
                hpp += ","
            comment = rxn.get('comment', '')
            if comment:
                hpp += f"  //  {comment}"
            hpp += "\n"

        hpp += "  };\n\n"

        # Equilibrium constants
        hpp += f"constexpr CArrayWrapper<double, {num_reactions}> equilibriumConstants = \n  {{ \n"
        for i, rxn in enumerate(self.reactions):
            val = self.format_scientific(rxn['equilibriumConstant'])
            hpp += f"    {val:>12}"
            if i < num_reactions - 1:
                hpp += ","
            comment = rxn.get('comment', '')
            if comment:
                hpp += f"   //  {comment}"
            hpp += "\n"
        hpp += "  };\n\n"

        # Forward rates
        hpp += f"constexpr CArrayWrapper<double, {num_reactions}> forwardRates = \n  {{ \n"
        for i, rxn in enumerate(self.reactions):
            val = self.format_scientific(rxn.get('forwardRate', 1.0e10))
            hpp += f"    {val:>12}"
            if i < num_reactions - 1:
                hpp += ","
            comment = rxn.get('comment', '')
            if comment:
                hpp += f"   //  {comment}"
            hpp += "\n"
        hpp += "  };\n\n"

        # Reverse rates
        hpp += f"constexpr CArrayWrapper<double, {num_reactions}> reverseRates = \n  {{ \n"
        for i, rxn in enumerate(self.reactions):
            val = self.format_scientific(rxn.get('reverseRate', 1.0e10))
            hpp += f"    {val:>12}"
            if i < num_reactions - 1:
                hpp += ","
            comment = rxn.get('comment', '')
            if comment:
                hpp += f"   //  {comment}"
            hpp += "\n"
        hpp += "  };\n\n"

        # Mobile species flags
        hpp += f"constexpr CArrayWrapper<int, {num_reactions}> mobileSpeciesFlag = \n  {{ \n"
        for i, rxn in enumerate(self.reactions):
            flag = rxn.get('mobileFlag', 1)
            hpp += f"    {flag}"
            if i < num_reactions - 1:
                hpp += ","
            comment = rxn.get('comment', '')
            if comment:
                hpp += f"   //  {comment}"
            hpp += "\n"
        hpp += "  };\n"

        hpp += f"""}}

  using {self.namespace}SystemType = reactionsSystems::MixedReactionsParameters< double, int, signed char, {num_species}, {num_reactions}, {num_equilibrium} >;

  constexpr {self.namespace}SystemType {self.namespace}System( {self.namespace}::stoichMatrix, {self.namespace}::equilibriumConstants, {self.namespace}::forwardRates, {self.namespace}::reverseRates, {self.namespace}::mobileSpeciesFlag );

// *****UNCRUSTIFY-ON******
}} // namespace geochemistry
}} // namespace hpcReact
"""

        return hpp


def scan_reaction_json_files(reaction_dir: str) -> List[str]:
    """Find all .json files in reactions/data/ directory"""
    data_dir = os.path.join(reaction_dir, 'data')
    if not os.path.exists(data_dir):
        return []

    return [
        os.path.join(data_dir, f)
        for f in os.listdir(data_dir)
        if f.endswith('.json')
    ]


def add_system_to_registry(registry_file: str, system_hpp: str, system_type: str):
    """Add a new system to existing GeochemicalSystems.hpp without regenerating"""

    if not os.path.exists(registry_file):
        print(f"Error: {registry_file} not found")
        return False

    with open(registry_file, 'r') as f:
        lines = f.readlines()

    # Check if system already exists
    content = ''.join(lines)
    if f'#include "{system_hpp}"' in content:
        print(f"System {system_hpp} already in registry")
        return False

    if system_type in content:
        print(f"System type {system_type} already in registry")
        return False

    # Find where to insert the include (before #include <variant>)
    include_insert_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith('#include') and '.hpp"' in line and '<variant>' not in line:
            include_insert_idx = i + 1
        elif line.strip() == '#include <variant>':
            break

    if include_insert_idx is None:
        print("Error: Could not find include section")
        return False

    # Insert the new include
    lines.insert(include_insert_idx, f'#include "{system_hpp}"\n')

    # Find the variant closing (look for the line with >; that ends the variant)
    variant_end_idx = None
    in_variant = False
    for i, line in enumerate(lines):
        if 'using systemTypes = std::variant<' in line:
            in_variant = True
        elif in_variant and '>;' in line:
            variant_end_idx = i
            break

    if variant_end_idx is None:
        print("Error: Could not find variant closing")
        return False

    # Find the last type line (line before >;)
    last_type_idx = variant_end_idx - 1
    while last_type_idx >= 0 and not lines[last_type_idx].strip():
        last_type_idx -= 1

    # Check if last type ends with comma, if not add it
    if not lines[last_type_idx].rstrip().endswith(','):
        lines[last_type_idx] = lines[last_type_idx].rstrip() + ',\n'

    # Insert new type with same indentation as previous line
    last_line = lines[last_type_idx]
    indent = len(last_line) - len(last_line.lstrip())
    indent_str = last_line[:indent]
    lines.insert(last_type_idx + 1, f'{indent_str}{system_type},\n')

    # Write back
    with open(registry_file, 'w') as f:
        f.writelines(lines)

    print(f"Added {system_type} to {registry_file}")
    return True


def generate_registry_header(reaction_dir: str, output_file: str):
    """Auto-generate GeochemicalSystems.hpp with all detected systems"""

    # Find all .hpp files (both manual and generated)
    geochemistry_dir = os.path.join(reaction_dir, 'geochemistry')
    if not os.path.exists(geochemistry_dir):
        print(f"Warning: {geochemistry_dir} not found")
        return

    hpp_files = [
        f for f in os.listdir(geochemistry_dir)
        if f.endswith('.hpp') and f not in ['GeochemicalSystems.hpp', 'Parameters.hpp']
    ]

    # Extract namespace names from JSON files or .hpp files
    systems = []
    for hpp_file in hpp_files:
        basename = os.path.splitext(hpp_file)[0]
        # Convert filename to namespace (e.g., "Ultramafics.hpp" -> "ultramafics")
        namespace = basename.lower()
        systems.append({
            'file': hpp_file,
            'basename': basename,
            'namespace': namespace,
            'typeName': f'{namespace}SystemType'
        })

    year = datetime.now().year
    registry = f"""/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: (BSD-3-Clause)
 *
 * Copyright (c) {year}- Lawrence Livermore National Security LLC
 * All rights reserved
 *
 * See top level LICENSE files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#pragma once

"""

    # Include all system headers
    for sys in systems:
        registry += f'#include "{sys["file"]}"\n'

    registry += """
#include <variant>

namespace hpcReact
{

namespace geochemistry
{
using systemTypes = std::variant< """

    # Add all system types to variant
    registry += ",\n                                  ".join(s['typeName'] for s in systems)

    registry += """ >;

} // namespace geochemistry
} // namespace hpcReact
"""

    with open(output_file, 'w') as f:
        f.write(registry)

    print(f"Generated {output_file} with {len(systems)} systems:")
    for sys in systems:
        print(f"  - {sys['basename']}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate constexpr C++ reaction system from JSON',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
JSON Input Format (Recommended - Separate Primary/Secondary Species):
{
  "systemName": "Carbonate System",
  "namespace": "carbonate",
  "secondarySpecies": ["OH-", "CO2", "CO3-2", "CaHCO3+"],
  "primarySpecies": ["H+", "HCO3-", "Ca+2"],
  "numEquilibriumReactions": 3,
  "reactions": [
    {
      "comment": "OH- + H+ = H2O",
      "stoichiometry": [-1, 0, 0, 0, -1, 0, 0],
      "equilibriumConstant": 9.89e13,
      "forwardRate": 1.0e10,
      "reverseRate": 1.0e10,
      "mobileFlag": 1,
      "isEquilibrium": true
    }
  ]
}

Notes:
- Species order in stoichiometry: [secondary..., primary...]
- isEquilibrium: true for fast equilibrium reactions, false for kinetic
- numEquilibriumReactions: number of equilibrium reactions (optional, auto-detected)
- mobileFlag: 1 for mobile species, 0 for immobile
- forwardRate and reverseRate: reaction rate constants

Build System Integration:
Add to CMakeLists.txt to auto-generate systems at build time.
        """
    )

    parser.add_argument('input', nargs='?', help='Input JSON file (or use --scan)')
    parser.add_argument('-o', '--output', help='Output .hpp file')
    parser.add_argument('--scan', metavar='DIR',
                       help='Scan reactions/data/ directory and generate all systems')
    parser.add_argument('--generate-registry', metavar='REACTION_DIR',
                       help='Generate GeochemicalSystems.hpp registry file')
    parser.add_argument('--add-to-registry', action='store_true',
                       help='Add new system to existing GeochemicalSystems.hpp instead of regenerating')

    args = parser.parse_args()

    if args.generate_registry:
        output = os.path.join(args.generate_registry, 'geochemistry', 'GeochemicalSystems.hpp')
        generate_registry_header(args.generate_registry, output)
        return

    if args.scan:
        json_files = scan_reaction_json_files(args.scan)
        if not json_files:
            print(f"No .json files found in {args.scan}/data/")
            return

        print(f"Found {len(json_files)} reaction JSON files:")
        for json_file in json_files:
            print(f"  Processing: {json_file}")
            with open(json_file, 'r') as f:
                config = json.load(f)

            gen = ReactionSystemGenerator(config)

            # Output to geochemistry/ directory
            basename = config['namespace'].capitalize()
            output_file = os.path.join(args.scan, 'geochemistry', f'{basename}.hpp')

            hpp_content = gen.generate_hpp()
            with open(output_file, 'w') as f:
                f.write(hpp_content)

            print(f"    -> Generated: {output_file}")
            print(f"       Species: {len(gen.species)}, Reactions: {len(gen.reactions)}")

        # Auto-generate registry
        print("\nGenerating registry header...")
        generate_registry_header(args.scan,
                                os.path.join(args.scan, 'geochemistry', 'GeochemicalSystems.hpp'))

    elif args.input:
        if not args.output:
            parser.error("--output required when processing single file")

        with open(args.input, 'r') as f:
            config = json.load(f)

        gen = ReactionSystemGenerator(config)
        hpp_content = gen.generate_hpp()

        with open(args.output, 'w') as f:
            f.write(hpp_content)

        print(f"Generated {args.output}")
        print(f"  Species: {len(gen.species)}")
        print(f"  Reactions: {len(gen.reactions)}")

        # If --add-to-registry flag is set, add to GeochemicalSystems.hpp
        if args.add_to_registry:
            # Determine registry file location
            # Assume output is in geochemistry/ directory
            output_dir = os.path.dirname(args.output)
            registry_file = os.path.join(output_dir, 'GeochemicalSystems.hpp')

            basename = os.path.basename(args.output)
            namespace = config['namespace']
            type_name = f'{namespace}SystemType'

            add_system_to_registry(registry_file, basename, type_name)

    else:
        parser.error("Either INPUT file or --scan DIR required")


if __name__ == '__main__':
    main()
