# Automated Reaction System Generation

## Overview

Define geochemical reaction systems in JSON instead of manually writing C++ matrices. Files are automatically generated during CMake configuration.

## Quick Start

**1. Create JSON file** (e.g., `my_system.json`):
```json
{
  "systemName": "My Carbonate System",
  "namespace": "mycarbonate",
  "secondarySpecies": ["OH-", "CO2"],
  "primarySpecies": ["H+", "HCO3-", "Ca+2"],
  "reactions": [
    {
      "equation": "OH- + H+ = H2O",
      "equilibriumConstant": 9.89e13,
      "forwardRate": 1.4e11,
      "reverseRate": 1.43e-3
    },
    {
      "equation": "CO2 + H2O = H+ + HCO3-",
      "equilibriumConstant": 4.42e-7,
      "forwardRate": 0.039,
      "reverseRate": 8.92e4
    },
    {
      "equation": "CaCO3 + H+ = Ca+2 + HCO3-",
      "equilibriumConstant": 5.16e1,
      "forwardRate": 1.55e-6,
      "reverseRate": 3.0e-8
    }
  ]
}
```

This example has 2 equilibrium reactions (first 2) + 1 kinetic reaction (mineral dissolution).

**2. Run configuration:**
```bash
python3 scripts/config-build.py -hc hostconfigs/your-config.cmake -bt Release
```

JSON is automatically validated, then `My_system.hpp` is generated.

**3. Build:**
```bash
cd build-<hostname>-<buildtype>  # e.g., build-macOC_arm-release
make
```

Done! Your reaction system is now compiled into the HPCReact library.

## JSON Format

### Required Fields
```json
{
  "systemName": "Human-readable name",
  "namespace": "lowercasecppname",
  "secondarySpecies": ["Species", "from", "equilibrium"],
  "primarySpecies": ["Basis", "species"],
  "reactions": [...]
}
```

### Reaction Format
```json
{
  "equation": "A + B = C + D",
  "equilibriumConstant": 1.23e10,
  "forwardRate": 1.0e10,
  "reverseRate": 1.0e-3
}
```

### Key Rules

- **Number of equilibrium reactions = Number of secondary species**
- First N reactions (N = # secondary) are equilibrium (fast)
- Remaining reactions are kinetic (time-dependent)
- Species order in equation: `[secondary..., primary...]`
- **Pure kinetic systems:** `secondarySpecies: []` → all reactions are kinetic

### Equation Syntax

- Separator: ` + ` (with spaces)
- Coefficients: `2Cl-` or `2 Cl-`
- Water: `H2O` (auto-skipped)
- Solids: Not in species lists (auto-skipped, e.g., `CaCO3` in kinetic reactions)
- Charges: Part of name (`H+`, `Ca+2`, `SO4-2`)

## Examples

- **[example_simple.json](example_simple.json)** - Minimal (1 equilibrium reaction)
- **[example_carbonate.json](example_carbonate.json)** - Mixed system (5 equilibrium + 1 kinetic)
- **[example_pure_kinetic.json](example_pure_kinetic.json)** - Pure kinetic (2 kinetic only)

## Tools

**Validate before configuring:**
```bash
python3 ../../../scripts/validate_reaction_json.py my_system.json
```

**Manual generation (for testing):**
```bash
python3 ../../../scripts/generate_reaction_system.py \
    my_system.json -o ../geochemistry/My_system.hpp
```

## File Locations

```
HPCReact/
├── scripts/
│   ├── generate_reaction_system.py   # Generator
│   └── validate_reaction_json.py     # Validator
└── src/reactions/
    ├── data/                          # JSON definitions (you edit)
    │   ├── example_simple.json
    │   └── example_carbonate.json
    └── geochemistry/                  # Generated .hpp (auto-created)
        ├── GeochemicalSystems.hpp     # Registry (auto-updated)
        ├── Example_carbonate.hpp      # Auto-generated
        └── Ultramafics.hpp            # Manual (existing)
```

## Workflow

```
Edit JSON → config-build.py → .hpp generated → make → done
```

Every time you:
- Add new JSON → re-run `config-build.py` → new `.hpp` created
- Edit existing JSON → re-run `config-build.py` → `.hpp` updated
- `GeochemicalSystems.hpp` registry always auto-updated

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "Unknown species 'X'" | Check spelling/charges: `Ca+2` not `Ca2+` |
| "Need at least N reactions" | Need ≥ N reactions for N secondary species |
| Changes not reflected | Re-run `config-build.py` (not just `make`) |