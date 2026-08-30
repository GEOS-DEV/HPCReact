# EQ3/6 reference data

References for the hardcoded values in the following tests that carry `EQ36` or `Bdot` in their name:

| Test | Quantity compared | Reference |
|---|---|---|
| `../testCarbonateActivityVsEQ36.cpp` | log10(gamma) of all 16 species | the log gamma column of the `carbonate.3o` species distribution table |
| `testcarbonateSystem_Bdot` in `../testGeochemicalEquilibriumReactions.cpp` | molality of the 7 primary species at equilibrium | the molality column of the `carbonate.3o` species distribution table |
| `computeReactionRatesVsEQ36_carbonateSystem_Bdot` in `../testGeochemicalKineticReactions.cpp` | the calcite dissolution rate `k*A*(1 - Q/K)` | the calcite saturation state `log Q/K` from `carbonate.3o` fed into the TST law |
| `testTimeStep_carbonateSystem_Bdot` in `../testGeochemicalMixedReactions.cpp` | molality of the 7 primary species after 10 s of calcite dissolution | the molality column of the `calcite.6o` species distribution table |

| File | What it is |
|---|---|
| `cmpHPCReact.d0` | EQ3/6 `data0.com.V8.R6`, with the 25 C Debye-Huckel A and B replaced by the values `HPCReact` derives from physical constants |
| `carbonate.3i` | EQ3NR input: the carbonate brine of `testcarbonateSystemAllEquilibrium`, with species outside the 17-species model suppressed |
| `carbonate.3o` | EQ3NR output |
| `calcite.6i` | EQ6 input: the `carbonate.3o` pickup, reacted with calcite under the TST rate law at a constant 100 cm2 (0.01 m2 in `HPCReact`) for 10 s |
| `calcite.6o` | EQ6 output |

## Regenerating

```bash
eqpt cmpHPCReact.d0          # writes cmpHPCReact.d1
eq3nr cmpHPCReact.d1 carbonate.3i
eq6   cmpHPCReact.d1 calcite.6i
```

`eqpt`, `eq3nr` and `eq6` are built from <https://github.com/39alpha/eq3_6>, which packages LLNL
EQ3/6 version 8.0a with a Make-based build.
