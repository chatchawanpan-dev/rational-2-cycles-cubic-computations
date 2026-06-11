#!/usr/bin/env sage
"""
SageMath verification script for the computations cited in the manuscript.

This script is intentionally limited to the SageMath computations used in the
updated manuscript:
  * the auxiliary rank-zero elliptic curve in the 3-torsion argument;
  * the duplication formulas for the distinguished section P_t;
  * the division-polynomial identities used in the torsion and rank arguments.

Run:
    sage sage_verification/verify_sage_identities.sage
"""

from sage.all import *


def header(title):
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def show_identity(name, computed, expected):
    print("")
    print(name)
    print("computed =", computed)
    print("expected =", expected)
    assert computed == expected
    print("status   = OK")


def verify_auxiliary_3torsion_curve():
    header("Auxiliary elliptic curve in the 3-torsion argument")

    E = EllipticCurve(QQ, [0, 0, 0, -39, -70])
    rank_bounds = E.rank_bounds()
    torsion = E.torsion_subgroup()
    torsion_points = sorted(E.torsion_points(), key=str)

    print("E': y^2 = x^3 - 39*x - 70")
    print("rank_bounds =", rank_bounds)
    print("rank =", E.rank())
    print("gens =", E.gens())
    print("torsion =", torsion)
    print("torsion_points =", torsion_points)

    assert rank_bounds == (0, 0)
    assert E.rank() == 0
    assert tuple(torsion.invariants()) == (2, 2)
    assert set(torsion_points) == {
        E(0),
        E((-5, 0)),
        E((-2, 0)),
        E((7, 0)),
    }
    print("status   = OK")


def verify_distinguished_section_formulas():
    header("Duplication formulas and division-polynomial checks")

    R.<t> = PolynomialRing(QQ)
    K = Frac(R)
    E = EllipticCurve(K, [0, 0, 0, -QQ(16) / 3, 16 * t^2 + QQ(128) / 27])
    P = E([K(QQ(4) / 3), -K(4) * t])

    P2 = 2 * P
    P4 = 4 * P
    expected_P2 = E([K(-QQ(8) / 3), K(4) * t])
    expected_P4 = E([
        K(QQ(16) / 3 * t^2 + 4) / t^2,
        K(-4 * t^4 - 16 * t^2 - 8) / t^3,
    ])

    print("E_t': y^2 = x^3 - 16/3*x + (16*t^2 + 128/27)")
    print("P =", P)
    print("2P =", P2)
    print("4P =", P4)
    assert P2 == expected_P2
    assert P4 == expected_P4
    print("duplication status = OK")

    psi3 = E.division_polynomial(3)
    psi5 = E.division_polynomial(5)
    psi7 = E.division_polynomial(7)
    psi9 = E.division_polynomial(9)

    xi1 = K(QQ(4) / 3)
    xi2 = K(-QQ(8) / 3)
    xi4 = K(4) / t^2 + K(QQ(16) / 3)

    show_identity("psi_3(4/3)", psi3(xi1), 256 * t^2)
    show_identity("psi_5(4/3)", psi5(xi1), -ZZ(2)^24 * t^6 * (t^2 + 1))
    show_identity("psi_7(4/3)", psi7(xi1), ZZ(2)^48 * t^12 * (t^4 - t^2 - 1))
    show_identity(
        "psi_9(4/3)",
        psi9(xi1),
        ZZ(2)^80 * t^20 * (3 * t^6 + 4 * t^4 + 3 * t^2 + 1),
    )
    show_identity("psi_3(-8/3)", psi3(xi2), -256 * (2 * t^2 + 1))
    show_identity(
        "psi_5(-8/3)",
        psi5(xi2),
        -ZZ(2)^24 * (t^2 + 1) * (t^6 - 5 * t^4 - 5 * t^2 - 1),
    )
    show_identity(
        "psi_3(4/t^2 + 16/3)",
        psi3(xi4),
        256 * (2 * t^2 + 1) * (2 * t^8 + 4 * t^6 + 10 * t^4 + 10 * t^2 + 3) / t^8,
    )


if __name__ == "__main__":
    verify_auxiliary_3torsion_curve()
    verify_distinguished_section_formulas()
    print("")
    print("All SageMath manuscript checks completed successfully.")
