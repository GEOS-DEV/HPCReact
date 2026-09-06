#!/usr/bin/env python3
"""
Independent reference values for the testEquilibriumReactions _Bdot cases.

The reaction-extent formulation parameterizes composition by how far each reaction has run,

    c_i( xi ) = c0_i + sum_r  s_ri * xi_r

which satisfies mole balance by construction, and then drives every reaction quotient onto its
equilibrium constant. Writing the residual as a single sum rather than a ratio of two products,

    residual_a = sum_i s_ai * ln( a_i )  -  ln( K_a )

this script evaluates it and its Jacobian at xi = 0, then runs its own Newton iteration to the
equilibrium composition. Expected values must not be copied from HPCReact's own output, or the
tests degenerate into regression checks that cannot detect a wrong activity model.

The B-dot machinery is imported from bdotKineticReference.py so the two scripts cannot drift.

Usage:  python3 src/reactions/exampleSystems/unitTests/bdotEquilibriumExtentReference.py
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bdotKineticReference import N, STOICH, CONC, activities  # noqa: E402

# Equilibrium constants for the two reactions (bulkGeneric::simpleMixedEquilibriumConstants).
K_EQ = [1.0, 1.0]
R = len(K_EQ)


def composition(xi):
    """c_i = c0_i + sum_r s_ri * xi_r."""
    return [CONC[i] + sum(STOICH[r][i] * xi[r] for r in range(R)) for i in range(N)]


def residual_and_jacobian(xi):
    conc = composition(xi)
    _, act, dact_dc = activities(conc)

    residual = [0.0] * R
    jacobian = [[0.0] * R for _ in range(R)]
    for a in range(R):
        residual[a] = sum(
            STOICH[a][i] * math.log(act[i]) for i in range(N)
        ) - math.log(K_EQ[a])

        for b in range(R):
            jacobian[a][b] = sum(
                STOICH[a][i] / act[i]
                * sum(dact_dc[i][j] * STOICH[b][j] for j in range(N))
                for i in range(N)
                if STOICH[a][i] != 0
            )

    return conc, residual, jacobian


def solve_2x2(matrix, rhs):
    det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    return [
        (rhs[0] * matrix[1][1] - matrix[0][1] * rhs[1]) / det,
        (matrix[0][0] * rhs[1] - rhs[0] * matrix[1][0]) / det,
    ]


def newton_to_equilibrium(max_iterations=200, tolerance=1.0e-14):
    """Newton on the extents, with the same positivity damping the C++ applies."""
    xi = [0.0] * R
    for iteration in range(max_iterations):
        _, residual, jacobian = residual_and_jacobian(xi)
        norm = math.sqrt(sum(v * v for v in residual))
        if norm < tolerance:
            return xi, iteration, norm

        step = solve_2x2(jacobian, [-v for v in residual])

        scale = 1.0
        for i in range(N):
            current = CONC[i] + sum(STOICH[r][i] * xi[r] for r in range(R))
            delta = sum(STOICH[r][i] * step[r] for r in range(R))
            if current + delta < 1.0e-30:
                damped = (1.0e-30 - current) / delta
                if damped < scale:
                    scale = 0.9 * damped

        xi = [xi[r] + scale * step[r] for r in range(R)]

    raise RuntimeError("Newton did not converge")


def main():
    _, residual, jacobian = residual_and_jacobian([0.0] * R)

    print("computeResidualAndJacobianTest_Bdot:")
    print()
    print("  double const expectedResiduals[] =")
    print("  { " + ", ".join(f"{v:.17g}" for v in residual) + " };")
    print("  double const expectedJacobian[2][2] =")
    print("  { { " + ", ".join(f"{v:.17g}" for v in jacobian[0]) + " },")
    print("    { " + ", ".join(f"{v:.17g}" for v in jacobian[1]) + " } };")
    print()
    asymmetry = jacobian[0][1] - jacobian[1][0]
    print(f"  Jacobian asymmetry J01 - J10 = {asymmetry:.12e}")
    print("  (zero only for the Identity model; a Cholesky solve is invalid when it is not)")
    print()

    xi, iterations, norm = newton_to_equilibrium()
    conc, residual, _ = residual_and_jacobian(xi)

    print(f"testEnforceEquilibrium_Bdot:  converged in {iterations} iterations, "
          f"||residual|| = {norm:.3e}")
    print("  Q/K at the solution = " + ", ".join(f"{math.exp(v):.15f}" for v in residual))
    print()
    print("  double const expectedSpeciesConcentrations[5] =")
    print("  { " + ", ".join(f"{v:.17g}" for v in conc) + " };")


if __name__ == "__main__":
    main()