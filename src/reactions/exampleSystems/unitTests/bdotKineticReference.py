#!/usr/bin/env python3
"""
Independent reference values for testKineticReactions.computeReactionRatesTest_..._Bdot.

The expected rates and Jacobian in that test must not be copied from HPCReact's own output, or the
test degenerates into a regression check that cannot detect a wrong activity model. This script
re-derives them from the B-dot equations directly, and cross-checks the analytic Jacobian against a
central difference of the rates.

System and parameters mirror bulkGeneric::simpleKineticTestRateParams and
bulkGeneric::simpleActivityTestParams in src/reactions/exampleSystems/BulkGeneric.hpp; the physical
constants mirror src/constitutive/activity/DebyeHuckel.hpp.

Usage:  python3 src/reactions/exampleSystems/unitTests/bdotKineticReference.py
"""

import math

# ---- physical constants (DebyeHuckel.hpp) ----
PI = 3.141592653589793e00
E0 = 8.854187812800001e-12
ECH = 1.602176634000000e-19
KB = 1.380649000000000e-23
NA = 6.022140760000000e23

# ---- water properties and scale factors (Bdot.hpp) ----
RHO_W = 997.0479  # kg/m3
EPS_R = 78.54  # dimensionless
T_K = 298.15  # K
LN10 = 2.302585092994046e00
M_PER_ANGSTROM = 1.0e-10

A_GAMMA_LN = (
    ECH**3
    * math.sqrt(2.0 * PI * NA * RHO_W)
    / (4.0 * PI * E0 * EPS_R * KB * T_K) ** 1.5
)
A = A_GAMMA_LN / LN10  # log10 scale, as Bdot.hpp uses it
B = math.sqrt(2.0 * NA * RHO_W * ECH**2 / (E0 * EPS_R * KB * T_K)) * M_PER_ANGSTROM

# ---- system definition (BulkGeneric.hpp) ----
N = 5
CHARGE = [2.0, -1.0, 1.0, 0.0, -1.0]
ION_SIZE = [4.0, 3.5, 3.5, 3.5, 3.5]  # Angstrom
BDOT = [0.1, 0.1, 0.1, 0.0, 0.1]

STOICH = [[-2, 1, 1, 0, 0],
          [0, 0, -1, -1, 2]]
K_FORWARD = [1.0, 0.5]
K_REVERSE = [1.0, 0.5]

CONC = [1.0, 1.0e-16, 0.5, 1.0, 1.0e-16]

DI_DC = [0.5 * z * z for z in CHARGE]


def activity_coefficients(conc):
    """ln(gamma_i) and d ln(gamma_i)/dc_j for the B-dot model."""
    ionic_strength = sum(conc[i] * DI_DC[i] for i in range(N))
    sqrt_i = math.sqrt(ionic_strength)

    ln_gamma = [0.0] * N
    dln_gamma_dc = [[0.0] * N for _ in range(N)]
    for i in range(N):
        denom = 1.0 + B * ION_SIZE[i] * sqrt_i
        log10_dh = -A * CHARGE[i] ** 2 * sqrt_i / denom
        dlog10_dh_di = -0.5 * A * CHARGE[i] ** 2 / (sqrt_i * denom * denom)

        ln_gamma[i] = (log10_dh + BDOT[i] * ionic_strength) * LN10
        dln_gamma_di = LN10 * (dlog10_dh_di + BDOT[i]) if ionic_strength > 0.0 else 0.0
        for j in range(N):
            dln_gamma_dc[i][j] = dln_gamma_di * DI_DC[j]

    return ionic_strength, ln_gamma, dln_gamma_dc


def activities(conc):
    """a_i = c_i * gamma_i, and d a_i / d c_j."""
    _, ln_gamma, dln_gamma_dc = activity_coefficients(conc)
    gamma = [math.exp(v) for v in ln_gamma]

    act = [conc[i] * gamma[i] for i in range(N)]
    dact_dc = [
        [
            (gamma[i] if i == j else 0.0) + conc[i] * gamma[i] * dln_gamma_dc[i][j]
            for j in range(N)
        ]
        for i in range(N)
    ]
    return gamma, act, dact_dc


def rates_and_jacobian(conc):
    """Elementary-law rates and d(rate_r)/d(c_j)."""
    gamma, act, dact_dc = activities(conc)

    rates = []
    drate_dc = []
    for r, stoich_row in enumerate(STOICH):
        forward = 1.0
        reverse = 1.0
        for i, s in enumerate(stoich_row):
            if s < 0:
                forward *= act[i] ** (-s)
            elif s > 0:
                reverse *= act[i] ** s
        rates.append(K_FORWARD[r] * forward - K_REVERSE[r] * reverse)

        drate_da = [0.0] * N
        for k, s in enumerate(stoich_row):
            if s < 0:
                drate_da[k] = K_FORWARD[r] * (-s) * forward / act[k]
            elif s > 0:
                drate_da[k] = -K_REVERSE[r] * s * reverse / act[k]

        drate_dc.append(
            [sum(drate_da[k] * dact_dc[k][j] for k in range(N)) for j in range(N)]
        )

    return gamma, rates, drate_dc


def rates_only(conc):
    return rates_and_jacobian(conc)[1]


def main():
    ionic_strength, _, _ = activity_coefficients(CONC)
    gamma, rates, drate_dc = rates_and_jacobian(CONC)

    print(f"A (log10 scale) = {A!r}")
    print(f"B               = {B!r}")
    print(f"ionic strength  = {ionic_strength!r}")
    print("gamma           = " + ", ".join(f"{v:.6f}" for v in gamma))
    print()
    print("  double const expectedReactionRates[] =")
    print("  { " + ", ".join(f"{v:.17g}" for v in rates) + " };")
    print()
    print("  double const expectedReactionRatesDerivatives[][5] =")
    for r in range(len(STOICH)):
        lead = "  { {" if r == 0 else "    {"
        tail = "}," if r == 0 else "} };"
        print(lead + " " + ", ".join(f"{v:.17g}" for v in drate_dc[r]) + " " + tail)

    print()
    print("central-difference cross-check of the analytic Jacobian:")
    worst = 0.0
    for j in range(N):
        step = max(abs(CONC[j]), 1.0) * 1.0e-7
        plus = list(CONC)
        plus[j] += step
        minus = list(CONC)
        minus[j] -= step
        rates_plus, rates_minus = rates_only(plus), rates_only(minus)
        for r in range(len(STOICH)):
            fd = (rates_plus[r] - rates_minus[r]) / (2.0 * step)
            worst = max(worst, abs(fd - drate_dc[r][j]) / max(abs(fd), 1.0e-12))
    print(f"  worst relative error = {worst:.2e}  (central-difference truncation)")


if __name__ == "__main__":
    main()