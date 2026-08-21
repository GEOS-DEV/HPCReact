# EQ3NR reference data

Provenance for the hardcoded values in `../testCarbonateActivityVsEQ36.cpp`.

| File | What it is |
|---|---|
| `cmpHPCReact.d0` | EQ3/6 `data0.com.V8.R6`, with the 25 C Debye-Huckel A and B replaced by the values this code derives from physical constants |
| `carbonate.3i` | EQ3NR input: the carbonate brine of `testcarbonateSystemAllEquilibrium`, with species outside this 17-species model suppressed |
| `carbonate.3o` | EQ3NR output. The molalities and log10(gamma) in the test are read from its species distribution table |

## Regenerating

```bash
eqpt cmpHPCReact.d0          # writes cmpHPCReact.d1
eq3nr cmpHPCReact.d1 carbonate.3i
```

See `RUNNING_EQ3NR.md` in the eq3_6 distribution for the input format and option switches.

## The database patch

Only two lines differ from stock `cmp.d0`:

```
debye huckel a (adh)     0.5114  ->  0.5084960
debye huckel b (bdh)     0.3288  ->  0.3281557
```

EQ3/6 tabulates A and B; this code computes them from physical constants with
`rho_w = 997.0479, eps_r = 78.54, T = 298.15`. Patching the database rather than the code isolates
the comparison to the activity model itself, so any disagreement is a real difference rather than a
difference in A and B. Without the patch, A alone differs by 0.57%, which puts roughly 1% into
gamma for divalent ions.

## Precision floor

EQ3NR truncates log10(gamma) to four decimals rather than rounding, so its reported values are
biased low in magnitude by up to 1e-4. That, not any model difference, sets the tolerance in the
test. EQ3/6 also fits A and B over its temperature grid rather than reading the 25 C entry
directly, so the effective values differ from the tabulated ones in the fifth decimal.

## Known departure

EQ3/6 gives CO2(aq) a Drummond (1981) salting-out coefficient; this code gives every neutral
species gamma = 1. That species is excluded from the comparison and asserted separately.