#!/usr/bin/env sage
"""
Supplementary SageMath scripts for:
"Rational 2-Periodic Points for Cubic Polynomials"

This file reproduces:
1) Division-polynomial evaluations used in Section 4 / Appendix A.
2) Table tab:NFI-samples (Normal Form I).
3) Table tab:NFII (Normal Form II).
4) Table tab:Ea-arith (non-normal family).

Run:
    sage cubic_2periodic_supplementary.sage
"""

from sage.all import *
import sys


def q(x):
    return str(QQ(x))


def pair_q(x, y):
    return f"({q(x)},{q(y)})"


def torsion_invariants(E):
    return tuple(E.torsion_subgroup().invariants())


def torsion_plain(inv):
    if len(inv) == 0:
        return "trivial"
    if len(inv) == 1:
        return f"Z/{inv[0]}Z"
    return " x ".join(f"Z/{n}Z" for n in inv)


def torsion_tex(inv):
    if len(inv) == 0:
        return "trivial"
    if len(inv) == 1:
        return f"$\\Z/{inv[0]}\\Z$"
    return "$" + "\\times".join(f"\\Z/{n}\\Z" for n in inv) + "$"


def exact_rank(E):
    lo, hi = E.rank_bounds()
    if lo != hi:
        raise RuntimeError(f"Rank bounds not sharp for {E}: bounds=({lo},{hi})")
    return lo


def print_header(title):
    bar = "=" * len(title)
    print("\n" + bar)
    print(title)
    print(bar)


def print_latex_rows(rows):
    for row in rows:
        print(" & ".join(row) + r" \\")


def show_identity(name, computed, expected):
    print(f"\n{name}")
    print(f"computed = {computed}")
    print(f"expected = {expected}")
    assert computed == expected
    print("status   = OK")


def show_direct_evaluation(poly_name, poly, xi_label, xi_value):
    val = poly(xi_value)
    print(f"\nCompute {poly_name} at {xi_label}")
    print(f"{poly_name}({xi_label}) = {val}")
    print(f"factor({poly_name}({xi_label})) = {factor(val)}")
    return val


def verify_division_polynomials():
    print_header("Division polynomial checks")
    R.<a> = PolynomialRing(QQ)
    K = Frac(R)
    E = EllipticCurve(K, [0, 0, 0, -QQ(16) / 3, 16 * a^2 + QQ(128) / 27])

    psi3 = E.division_polynomial(3)
    psi5 = E.division_polynomial(5)
    psi7 = E.division_polynomial(7)
    psi9 = E.division_polynomial(9)

    print("Curve E_a' over Q(a): y^2 = x^3 - 16/3*x + (16*a^2 + 128/27)")
    print("\nExplicit division polynomials (odd indices):")
    print(f"psi_3(x) = {psi3}")
    print(f"psi_5(x) = {psi5}")
    print(f"psi_7(x) = {psi7}")
    print(f"psi_9(x) = {psi9}")

    xi1 = QQ(4) / 3
    lhs1 = show_direct_evaluation("psi_3", psi3, "4/3", xi1)
    rhs1 = 256 * a^2
    lhs2 = show_direct_evaluation("psi_5", psi5, "4/3", xi1)
    rhs2 = -ZZ(2)^24 * a^6 * (a^2 + 1)
    lhs3 = show_direct_evaluation("psi_7", psi7, "4/3", xi1)
    rhs3 = ZZ(2)^48 * a^12 * (a^4 - a^2 - 1)
    lhs4 = show_direct_evaluation("psi_9", psi9, "4/3", xi1)
    rhs4 = ZZ(2)^80 * a^20 * (3 * a^6 + 4 * a^4 + 3 * a^2 + 1)

    xi2 = -QQ(8) / 3
    lhs5 = show_direct_evaluation("psi_3", psi3, "-8/3", xi2)
    rhs5 = -256 * (2 * a^2 + 1)
    lhs6 = show_direct_evaluation("psi_5", psi5, "-8/3", xi2)
    rhs6 = -ZZ(2)^24 * (a^2 + 1) * (a^6 - 5 * a^4 - 5 * a^2 - 1)

    x4 = QQ(4) / a^2 + QQ(16) / 3
    lhs7 = show_direct_evaluation("psi_3", psi3, "4/a^2 + 16/3", x4)
    rhs7 = 256 * (2 * a^2 + 1) * (2 * a^8 + 4 * a^6 + 10 * a^4 + 10 * a^2 + 3) / a^8

    print("\nEvaluations used in the manuscript:")
    show_identity("At xi = 4/3: psi_3(4/3)", lhs1, rhs1)
    show_identity("At xi = 4/3: psi_5(4/3)", lhs2, rhs2)
    show_identity("At xi = 4/3: psi_7(4/3)", lhs3, rhs3)
    show_identity("At xi = 4/3: psi_9(4/3)", lhs4, rhs4)
    show_identity("At xi = -8/3: psi_3(-8/3)", lhs5, rhs5)
    show_identity("At xi = -8/3: psi_5(-8/3)", lhs6, rhs6)
    show_identity("At xi = 4/a^2 + 16/3: psi_3(4/a^2 + 16/3)", lhs7, rhs7)

def reproduce_table_nfi():
    print_header("Table tab:NFI-samples")

    rows_out = []

    symmetric_specs = [
        (QQ(2), QQ(1)),
        (QQ(1), QQ(2)),
        (QQ(4), QQ(1) / 2),
        (QQ(2), QQ(3)),
    ]
    for r, a in symmetric_specs:
        b = -1 - a * r^2 / 4
        x1 = r / 2
        x2 = -r / 2
        phi = lambda z: a * z^3 + b * z
        assert phi(x1) == x2 and phi(x2) == x1 and x1 != x2
        rows_out.append([
            "Symmetric",
            f"$r={q(r)}$",
            f"${q(a)}$",
            f"${q(b)}$",
            f"$({q(x1)},{q(x2)})$",
            f"$({q(phi(x1))},{q(phi(x2))})$",
        ])

    generic_specs = [
        (QQ(1), QQ(2)),
        (QQ(1), QQ(3)),
        (QQ(2), QQ(5)),
        (QQ(3) / 2, QQ(5) / 2),
    ]
    for u, v in generic_specs:
        assert u != 0 and v != 0 and v^2 != u^2
        a = 4 / (v^2 - u^2)
        b = -2 * (u^2 + v^2) / (v^2 - u^2)
        x1 = (u + v) / 2
        x2 = (u - v) / 2
        phi = lambda z: a * z^3 + b * z
        assert phi(x1) == x2 and phi(x2) == x1 and x1 != x2
        rows_out.append([
            "Generic",
            f"$(u,v)=({q(u)},{q(v)})$",
            f"${q(a)}$",
            f"${q(b)}$",
            f"$({q(x1)},{q(x2)})$",
            f"$({q(phi(x1))},{q(phi(x2))})$",
        ])

    print_latex_rows(rows_out)
    return rows_out


def reproduce_table_nfii():
    print_header("Table tab:NFII")

    # Exact row order from the manuscript.
    As = [
        QQ(-1), QQ(-1) / 2,
        QQ(-8), QQ(16), QQ(-48), QQ(72), QQ(1), QQ(-2), QQ(5), QQ(12),
        QQ(6), QQ(7), QQ(-9), QQ(14), QQ(21), QQ(25),
        QQ(-30), QQ(40), QQ(-47), QQ(53),
    ]

    expected_rank = {
        QQ(-1): 0, QQ(-1) / 2: 0,
        QQ(-8): 1, QQ(16): 1, QQ(-48): 1, QQ(72): 1, QQ(1): 1, QQ(-2): 1, QQ(5): 1, QQ(12): 1,
        QQ(6): 2, QQ(7): 2, QQ(-9): 2, QQ(14): 2, QQ(21): 2, QQ(25): 2,
        QQ(-30): 3, QQ(40): 3, QQ(-47): 3, QQ(53): 3,
    }

    expected_torsion = {
        QQ(-1): (5,),
        QQ(-1) / 2: (6,),
        QQ(-8): (2,),
        QQ(16): (2,),
        QQ(-48): (2,),
        QQ(72): (2,),
    }

    rows_out = []
    rank0_count = 0

    for A in As:
        assert A != 0 and 27 * A + 16 != 0
        E = EllipticCurve([0, 4 * A, 0, 0, 16 * A^4])  # model for E^(A)
        rank = exact_rank(E)
        inv = torsion_invariants(E)
        if rank == 0:
            rank0_count += 1

        exp_rank = expected_rank[A]
        exp_tor = expected_torsion.get(A, ())
        assert rank == exp_rank, (A, rank, exp_rank)
        assert inv == exp_tor, (A, inv, exp_tor)

        # Canonical rational point on E^(A) in this model:
        # T_A = (-4A,-4A^2), corresponding to (X,Y)=(1,-A).
        T = E([-4 * A, -4 * A^2])
        k = 1 if rank == 0 else -2
        Q = k * T
        x, y = Q.xy()

        # Recovery formulas from the manuscript in this model:
        # B = -x/(4A) - 2 - 16A^3/x^2
        # z1 = (-4A^2 - y)/(2Ax), z2 = (-4A^2 + y)/(2Ax)
        B = -x / (4 * A) - 2 - 16 * A^3 / x^2
        z1 = (-4 * A^2 - y) / (2 * A * x)
        z2 = (-4 * A^2 + y) / (2 * A * x)

        psi = lambda z: A * z^3 + B * z + 1
        assert psi(z1) == z2 and psi(z2) == z1 and z1 != z2

        rows_out.append([
            f"${q(A)}$",
            torsion_tex(inv),
            f"${rank}$",
            f"${q(B)}$",
            f"$({q(z1)},{q(z2)})$",
        ])

    print_latex_rows(rows_out)
    pct_rank0 = 100.0 * rank0_count / len(rows_out)
    print(f"\nrank-0 rows: {rank0_count}/{len(rows_out)} = {pct_rank0:.1f}%")
    return rows_out


def reproduce_table_non_normal():
    print_header("Table tab:Ea-arith")

    # Exact row order from the manuscript.
    a_values = [
        QQ(1) / 2, QQ(1), QQ(3) / 2, QQ(2), QQ(5) / 2, QQ(3), QQ(4), QQ(6), QQ(20),
        QQ(7) / 2, QQ(5), QQ(13) / 2, QQ(9), QQ(11), QQ(14), QQ(21) / 2,
        QQ(17) / 2, QQ(19), QQ(29), QQ(41) / 2,
    ]

    expected_rank = {
        QQ(1) / 2: 1, QQ(1): 1, QQ(3) / 2: 1, QQ(2): 1, QQ(5) / 2: 1, QQ(3): 1, QQ(4): 1, QQ(6): 1, QQ(20): 1,
        QQ(7) / 2: 2, QQ(5): 2, QQ(13) / 2: 2, QQ(9): 2, QQ(11): 2, QQ(14): 2, QQ(21) / 2: 2,
        QQ(17) / 2: 3, QQ(19): 3, QQ(29): 3, QQ(41) / 2: 3,
    }

    expected_torsion = {
        QQ(4): (2,),
        QQ(20): (2,),
    }

    rows_out = []
    nontrivial_torsion_count = 0

    for a in a_values:
        assert a != 0
        E = EllipticCurve([0, 4, 0, 0, 16 * a^2])  # E_a: y^2 = x^3 + 4x^2 + 16a^2
        rank = exact_rank(E)
        inv = torsion_invariants(E)
        if len(inv) > 0:
            nontrivial_torsion_count += 1

        exp_rank = expected_rank[a]
        exp_tor = expected_torsion.get(a, ())
        assert rank == exp_rank, (a, rank, exp_rank)
        assert inv == exp_tor, (a, inv, exp_tor)

        # As described in the manuscript: Q = [3]P_a with P_a = (0,-4a).
        P = E([0, -4 * a])
        Q = 3 * P
        X, Y = Q.xy()

        b = -X / 4 - 2 - 16 * a^2 / X^2
        x_plus = (-4 * a - Y) / (2 * X)
        x_minus = (-4 * a + Y) / (2 * X)

        f = lambda x: x^3 + b * x + a
        assert f(x_plus) == x_minus and f(x_minus) == x_plus and x_plus != x_minus

        rows_out.append([
            f"${q(a)}$",
            f"${q(b)}$",
            f"$({q(x_plus)},{q(x_minus)})$",
            torsion_tex(inv),
            f"${rank}$",
        ])

    print_latex_rows(rows_out)
    pct_nontrivial = 100.0 * nontrivial_torsion_count / len(rows_out)
    print(
        f"\nnontrivial torsion rows: {nontrivial_torsion_count}/{len(rows_out)}"
        f" = {pct_nontrivial:.1f}%"
    )
    return rows_out


if __name__ == "__main__":
    division_only = "--division-only" in sys.argv
    if division_only:
        verify_division_polynomials()
    else:
        verify_division_polynomials()
        reproduce_table_nfi()
        reproduce_table_nfii()
        reproduce_table_non_normal()
