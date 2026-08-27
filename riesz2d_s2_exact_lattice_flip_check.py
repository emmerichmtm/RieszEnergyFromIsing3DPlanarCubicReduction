#!/usr/bin/env python3
"""Lattice-sum constants, objective coefficient, and component-flip loads.

Companion to arXiv:2608.23506 (Riesz 2-energy subset selection in the plane).
Written during the independent soundness audit performed with Claude Fable 5
(Anthropic), 27 August 2026.  Uses only the Python standard library.

Checks performed (assertions; non-zero exit on failure):

  1. zeta constants: zeta(5)-1 < 0.037, Q6 <= zeta(3)^2/8 < 0.181, P6 < 0.020,
     and the U-turn sum U6 = sum_{k>=1} (k+1)(k^2+1)^-3 < 0.281 (eq. (46) of v3).
  2. The "elementary certification" arithmetic: the v2 bound
     sum_{k=2..5} k^-5 + int_5^inf x^-5 = 0.03706 does NOT certify 0.037;
     the v3 bound with k up to 6 does.
  3. Lemma 8.1: 19 b^4 < lambda < 37 b^4 for D = 2, 10, 100, 1000 (leading
     formula), and the exact rational lambda/b^4 = 29.270918... at D = 2,
     b = 1/500 reported in Appendix A.
  4. Lemma 9.1: absolute cross-cut coupling load for every cut type -- straight
     run, root junction (both branch configurations), and every cut inside the
     objective-box template including the unit-step U-turn -- stays below the
     paper's budget 95 b^4 (same tree) + 37 b^4 (cross tree), against a
     favorable edge of 160 b^4.  Both the leading-order |J| and the crude
     225 r^-6 bound used in the paper are reported.
"""
import math
from fractions import Fraction as Fr

from riesz2d_s2_exact_multipole_check import fourier  # noqa: E402


def zeta(s, N=200000):
    return sum(k ** -s for k in range(1, N)) + N ** (1 - s) / (s - 1)


def Jlead(R):
    """Leading spin-spin coefficient in units of b^4 for integer R."""
    r2 = R[0] ** 2 + R[1] ** 2
    s2 = 2 * R[0] * R[1] / r2
    return (384 * s2 ** 2 - 160) / r2 ** 3


def load(C, rest, favorable_edge):
    """Sum of |J| (leading) and of 225 r^-6 over cross-cut pairs C x rest."""
    tot = crude = 0.0
    for a in C:
        for c in rest:
            if {a, c} == set(favorable_edge):
                continue
            R = (c[0] - a[0], c[1] - a[1])
            tot += abs(Jlead(R))
            crude += 225 / (R[0] ** 2 + R[1] ** 2) ** 3
    return tot, crude


def main():
    print("== 1. Lattice-sum constants ==")
    z5m1 = zeta(5) - 1
    Q6 = sum((i * i + j * j) ** -3 for i in range(1, 3000) for j in range(1, 3000))
    P6 = zeta(5) - zeta(6)
    U6 = sum((k + 1) * (k * k + 1) ** -3 for k in range(1, 200000))
    print(f"  zeta(5)-1 = {z5m1:.6f}   Q6 = {Q6:.6f} (bound {zeta(3) ** 2 / 8:.6f})"
          f"   P6 = {P6:.6f}   U6 = {U6:.6f}")
    assert z5m1 < 0.037 and Q6 < 0.181 and P6 < 0.020 and U6 < 0.281
    budget = 225 * (0.037 + 2 * 0.181 + 0.020 + 0.001)
    uturn = 225 * (0.037 + 0.281 + 0.001)
    print(f"  (45) budget 225*(0.037+2*0.181+0.020+0.001) = {budget:.3f} < 95")
    print(f"  U-turn budget 225*(0.037+0.281+0.001) = {uturn:.3f} < 72")
    assert budget < 95 and uturn < 72

    print("\n== 2. Elementary certification arithmetic ==")
    v2 = sum(k ** -5 for k in range(2, 6)) + 5 ** -4 / 4
    v3 = sum(k ** -5 for k in range(2, 7)) + 6 ** -4 / 4
    u_partial = sum((k + 1) * (k * k + 1) ** -3 for k in range(1, 6))
    u_tail = (6 ** -5 + 6 ** -4 / 4) + (6 ** -6 + 6 ** -5 / 5)
    print(f"  v2 (k<=5, int from 5): {v2:.6f}  -> certifies 0.037? {v2 < 0.037}")
    print(f"  v3 (k<=6, int from 6): {v3:.6f}  -> certifies 0.037? {v3 < 0.037}")
    print(f"  U6 partial (k<=5) = {u_partial:.6f} < 0.2794; tail bound = {u_tail:.6f} < 0.0004")
    assert not (v2 < 0.037) and v3 < 0.037
    assert u_partial < 0.2794 and u_tail < 0.0004

    print("\n== 3. Objective coefficient lambda (Lemma 8.1, Appendix A) ==")
    for D in (2, 10, 100, 1000):
        lam = sum(Jlead((1 + i + j, 1)) for i in range(D + 1) for j in range(D + 1))
        halo = sum(abs(Jlead((1 + i + j, 1))) for i in range(D + 1) for j in range(D + 1)
                   if (i, j) != (0, 0))
        print(f"  D={D:5d}: lambda/b^4 = {lam:.6f}   |halo|/b^4 = {halo:.4f}")
        assert 19 < lam < 37 and halo < 8.4
    b = Fr(1, 500)
    lam_exact = sum(fourier((Fr(1 + i + j), Fr(1)), b)[3] for i in range(3) for j in range(3))
    print(f"  exact lambda/b^4 at D=2, b=1/500: {float(lam_exact / b ** 4):.6f}  (Appendix A: 29.270918)")
    assert abs(float(lam_exact / b ** 4) - 29.270918) < 1e-6

    print("\n== 4. Component-flip cross-cut loads (units of b^4) ==")
    D, N = 40, 120
    q = [(1 + j, 1) for j in range(D + 1)]          # T_v prefix; q_0 is the leaf end
    step = (1 + D, 0)
    ray = [(1 + D + m, 0) for m in range(1, N)]     # T_v continuation
    p = [(-i, 0) for i in range(N)]                 # T_u prefix and continuation
    cuts = {
        "unit step  (q_D | step)": (q, [step] + ray, (q[-1], step)),
        "q_{D-1} | q_D": (q[:-1], [q[-1], step] + ray, (q[-2], q[-1])),
        "q_{D-2} | q_{D-1}": (q[:-2], q[-2:] + [step] + ray, (q[-3], q[-2])),
        "q_0 | q_1": (q[:1], q[1:] + [step] + ray, (q[0], q[1])),
        "step | ray_1": (q + [step], ray, (step, ray[0])),
        "ray_1 | ray_2": (q + [step] + ray[:1], ray[1:], (ray[0], ray[1])),
    }
    worst = 0.0
    for name, (C, rest, e) in cuts.items():
        same, same_crude = load(C, rest, e)
        xt, xt_crude = load(C, p, (None, None))
        fav = abs(Jlead((e[1][0] - e[0][0], e[1][1] - e[0][1])))
        print(f"  {name:26s} same-tree {same:6.2f} (crude {same_crude:6.2f})"
              f"  cross-tree {xt:6.2f} (crude {xt_crude:6.2f})  favorable {fav:.1f}")
        assert same_crude < 95 and xt_crude < 37 and fav > 159
        worst = max(worst, same_crude + xt_crude)
    for name, others in [("root: right/up/left", [(0, j) for j in range(1, N)] + [(-j, 0) for j in range(1, N)]),
                         ("root: right/up/down", [(0, j) for j in range(1, N)] + [(0, -j) for j in range(1, N)])]:
        C = [(i, 0) for i in range(1, N)]
        same, crude = load(C, [(0, 0)] + others, ((1, 0), (0, 0)))
        print(f"  {name:26s} same-tree {same:6.2f} (crude {crude:6.2f})")
        assert crude < 95
        worst = max(worst, crude + 37)
    print(f"  worst total crude adverse load {worst:.2f} b^4 < 132.75 b^4 budget; favorable edge > 158.6 b^4")
    assert worst < 132.75

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
