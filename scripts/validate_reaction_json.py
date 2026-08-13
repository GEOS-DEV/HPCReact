#!/usr/bin/env python3
"""
Validate reaction system JSON files before generating C++ code
"""

import argparse
import json
import sys
from typing import Dict, List


class ReactionValidator:
    def __init__(self, config: Dict):
        self.config = config
        self.errors = []
        self.warnings = []

    def validate(self) -> bool:
        """Run all validation checks"""
        self._validate_required_fields()
        self._validate_species()
        self._validate_reactions()
        self._validate_namespace()

        return len(self.errors) == 0

    def _validate_required_fields(self):
        """Check all required fields are present"""
        required = ['systemName', 'namespace', 'primarySpecies', 'secondarySpecies', 'reactions']
        for field in required:
            if field not in self.config:
                self.errors.append(f"Missing required field: '{field}'")

    def _validate_species(self):
        """Validate species lists"""
        if 'primarySpecies' not in self.config or 'secondarySpecies' not in self.config:
            return

        primary = self.config['primarySpecies']
        secondary = self.config['secondarySpecies']

        if not isinstance(primary, list) or len(primary) == 0:
            self.errors.append("primarySpecies must be a non-empty list")

        if not isinstance(secondary, list):
            self.errors.append("secondarySpecies must be a list (can be empty for pure kinetic systems)")

        # If no secondary species, this is a pure kinetic system
        if len(secondary) == 0:
            self.warnings.append("Pure kinetic system detected (no secondary species = no equilibrium reactions)")

        # Check for duplicates
        all_species = primary + secondary
        if len(all_species) != len(set(all_species)):
            duplicates = [s for s in all_species if all_species.count(s) > 1]
            self.errors.append(f"Duplicate species: {', '.join(set(duplicates))}")

        # Check for common naming issues
        for sp in all_species:
            if ' ' in sp:
                self.warnings.append(f"Species '{sp}' contains space - this may cause parsing issues")

    def _validate_reactions(self):
        """Validate reaction definitions"""
        if 'reactions' not in self.config:
            return

        reactions = self.config['reactions']
        if not isinstance(reactions, list) or len(reactions) == 0:
            self.errors.append("reactions must be a non-empty list")
            return

        num_secondary = len(self.config.get('secondarySpecies', []))

        if num_secondary > 0 and len(reactions) < num_secondary:
            self.errors.append(
                f"Need at least {num_secondary} reactions (one per secondary species), "
                f"got {len(reactions)}"
            )

        # Validate each reaction
        for i, rxn in enumerate(reactions):
            self._validate_reaction(i, rxn, num_secondary)

    def _validate_reaction(self, idx: int, rxn: Dict, num_secondary: int):
        """Validate a single reaction"""
        rxn_label = f"Reaction {idx+1}"

        # Check required fields
        if 'equation' not in rxn:
            self.errors.append(f"{rxn_label}: Missing 'equation' field")
            return

        if 'equilibriumConstant' not in rxn:
            self.errors.append(f"{rxn_label}: Missing 'equilibriumConstant' field")

        # Check equation format
        equation = rxn['equation']
        if '=' not in equation:
            self.errors.append(f"{rxn_label}: Equation must contain '=' : {equation}")

        # Parse and validate species
        try:
            left, right = equation.split('=')
            all_species = self.config.get('primarySpecies', []) + self.config.get('secondarySpecies', [])

            for side in [left, right]:
                # Split by + surrounded by spaces to avoid splitting charges
                import re
                terms = re.split(r'\s+\+\s+', side)
                for term in terms:
                    term = term.strip()
                    if not term or term in ['H2O', 'H₂O']:
                        continue

                    # Remove leading coefficients (digits before capital letter)
                    match = re.match(r'^\d+\s*([A-Z].*)$', term)
                    if match:
                        term = match.group(1).strip()

                    # Check if species is known or likely a solid
                    if term not in all_species:
                        # Likely a solid/mineral - this is OK
                        if idx >= num_secondary:
                            # Kinetic reaction - solids expected
                            pass
                        else:
                            self.warnings.append(
                                f"{rxn_label}: Species '{term}' not in primary/secondary lists "
                                f"(OK if it's a solid/mineral)"
                            )
        except Exception as e:
            self.errors.append(f"{rxn_label}: Failed to parse equation: {e}")

        # Check rate constants
        forward = rxn.get('forwardRate', 1.0e10)
        reverse = rxn.get('reverseRate', 1.0e10)

        if forward <= 0:
            self.errors.append(f"{rxn_label}: forwardRate must be positive")
        if reverse <= 0:
            self.errors.append(f"{rxn_label}: reverseRate must be positive")

        # Warn about equilibrium/kinetic classification
        is_equilibrium = (idx < num_secondary)
        if is_equilibrium:
            # Equilibrium reactions should be fast
            if forward < 1e6 or reverse < 1e6:
                self.warnings.append(
                    f"{rxn_label}: Equilibrium reaction but has slow rates "
                    f"(forward={forward:.2e}, reverse={reverse:.2e})"
                )
        else:
            # Kinetic reactions should be slower
            if forward >= 1e9 and reverse >= 1e9:
                self.warnings.append(
                    f"{rxn_label}: Kinetic reaction but has very fast rates "
                    f"(forward={forward:.2e}, reverse={reverse:.2e}). "
                    f"Consider moving to equilibrium reactions."
                )

    def _validate_namespace(self):
        """Validate C++ namespace"""
        if 'namespace' not in self.config:
            return

        namespace = self.config['namespace']

        # Must be valid C++ identifier
        if not namespace.replace('_', '').isalnum():
            self.errors.append(f"namespace '{namespace}' is not a valid C++ identifier")

        if namespace[0].isdigit():
            self.errors.append(f"namespace '{namespace}' cannot start with a digit")

        if not namespace.islower():
            self.warnings.append(f"namespace '{namespace}' should be lowercase")

    def print_results(self):
        """Print validation results"""
        if self.errors:
            print("❌ ERRORS:")
            for error in self.errors:
                print(f"  - {error}")
            print()

        if self.warnings:
            print("⚠️  WARNINGS:")
            for warning in self.warnings:
                print(f"  - {warning}")
            print()

        if not self.errors and not self.warnings:
            print("✅ Validation passed! No errors or warnings.")
        elif not self.errors:
            print("✅ Validation passed with warnings.")
        else:
            print("❌ Validation failed.")


def main():
    parser = argparse.ArgumentParser(
        description='Validate reaction system JSON files'
    )
    parser.add_argument('json_file', help='JSON file to validate')
    parser.add_argument('--strict', action='store_true',
                       help='Treat warnings as errors')

    args = parser.parse_args()

    try:
        with open(args.json_file, 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        print(f"❌ Error: File not found: {args.json_file}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"❌ Error: Invalid JSON: {e}")
        sys.exit(1)

    validator = ReactionValidator(config)
    is_valid = validator.validate()
    validator.print_results()

    if not is_valid:
        sys.exit(1)

    if args.strict and validator.warnings:
        sys.exit(1)

    sys.exit(0)


if __name__ == '__main__':
    main()
