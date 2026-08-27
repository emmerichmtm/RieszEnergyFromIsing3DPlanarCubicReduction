#!/usr/bin/env python3
"""Exact-rational verification of the angular multipole law (Lemma 5.1).

Companion to arXiv:2608.23506 (Riesz 2-energy subset selection in the plane).
Written during the independent soundness audit performed with Claude Fable 5
(Anthropic), 27 August 2026.  Uses only the Python standard library.

At exponent s = 2 every pair energy of rational points is rational, so the
exact two-selector interaction F_R(sigma, tau) of eq. (14) and its Fourier
coefficients (16) can be evaluated with `fractions.Fraction` and compared
against the closed-form leading terms of Lemma 5.1 without any rounding.

Checks performed (each is an assertion; the script exits non-zero on failure):

  1. Figure 3 values at b = 0.12 for R = (0,1) and R = (1,1).
  2. J_R / (b^4 r^-6)  ->  384 sin^2(2 theta) - 160   as b -> 0,
     h_R / (b^2 r^-4)  ->  16 sin(2 theta),
     for the center vectors that occur in the construction.
  3. Exact structural identities: h^(1) = h^(2); F_R even in b (all odd orders
     vanish); one-spin field identically zero for axis-aligned R.
  4. The remainder |R_J| + |R_h| is O(t^6 r^-2) with coefficient < 4000,
     hence far below the uniform bound 1e6 t^5 r^-2 of eq. (19) for t <= 1e-8.
  5. Internal field coefficient B(eta) of eq. (23): beta = -b^-3/8 and the
     quadratic residual constant of eq. (32) is about 0.19 < 1.
  6. Sensitivity of Fourier coefficients to a diagonal-radius perturbation is
     far below the constant 100 of eq. (33).
"""
from fractions import Fraction as Fr

U_PLUS = (Fr(1), Fr(1))
U_MINUS = (Fr(1), Fr(-1))


def u(s):
    return U_PLUS if s == 1 else U_MINUS


def F(R, b, sigma, tau, eta1=Fr(0), eta2=Fr(0)):
    """Exact four-pair interaction (14); the + diagonal has radius b + eta."""
    total = Fr(0)
    r1 = b + eta1 if sigma == 1 else b
    r2 = b + eta2 if tau == 1 else b
    for e in (1, -1):
        for n in (1, -1):
            dx = R[0] + n * r2 * u(tau)[0] - e * r1 * u(sigma)[0]
            dy = R[1] + n * r2 * u(tau)[1] - e * r1 * u(sigma)[1]
            total += 1 / (dx * dx + dy * dy)
    return total


def fourier(R, b, eta1=Fr(0), eta2=Fr(0)):
    """Exact Fourier coefficients (A, h1, h2, J) of eq. (16)."""
    pp = F(R, b, 1, 1, eta1, eta2)
    pm = F(R, b, 1, -1, eta1, eta2)
    mp = F(R, b, -1, 1, eta1, eta2)
    mm = F(R, b, -1, -1, eta1, eta2)
    A = (pp + pm + mp + mm) / 4
    h1 = (pp + pm - mp - mm) / 4
    h2 = (pp - pm + mp - mm) / 4
    J = (pp - pm - mp + mm) / 4
    return A, h1, h2, J


def leading_J(R):
    """(384 sin^2(2 theta) - 160) / r^6, exactly rational for integer R."""
    r2 = R[0] ** 2 + R[1] ** 2
    s2 = 2 * R[0] * R[1] / r2
    return (384 * s2 ** 2 - 160) / r2 ** 3


def leading_h(R, b):
    """16 sin(2 theta) b^2 r^-4 + 128 sin(2 theta) b^4 r^-6."""
    r2 = R[0] ** 2 + R[1] ** 2
    s2 = 2 * R[0] * R[1] / r2
    return 16 * s2 * b ** 2 / r2 ** 2 + 128 * s2 * b ** 4 / r2 ** 3


def main():
    print("== 1. Figure 3 (b = 0.12) ==")
    b = Fr(12, 100)
    fig3 = {
        (0, 1): dict(J=-0.035776, Fpp=4.201188, Fpm=4.272740),
        (1, 1): dict(J=0.006145, Fpp=2.190833, Fpm=2.055895),
    }
    for Rint, ref in fig3.items():
        R = (Fr(Rint[0]), Fr(Rint[1]))
        J = fourier(R, b)[3]
        Fpp, Fpm = F(R, b, 1, 1), F(R, b, 1, -1)
        print(f"  R={Rint}: J={float(J):.6f}  F(+,+)={float(Fpp):.6f}  F(+,-)={float(Fpm):.6f}")
        assert abs(float(J) - ref["J"]) < 1e-6
        assert abs(float(Fpp) - ref["Fpp"]) < 1e-6
        assert abs(float(Fpm) - ref["Fpm"]) < 1e-6

    print("\n== 2-4. Leading coefficients, exact identities, remainder ==")
    vectors = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (3, 1), (5, 2), (7, 3)]
    max_rem6 = Fr(0)
    for Rint in vectors:
        R = (Fr(Rint[0]), Fr(Rint[1]))
        r2 = R[0] ** 2 + R[1] ** 2
        s2 = 2 * R[0] * R[1] / r2
        for b in (Fr(1, 1000), Fr(1, 10000)):
            A, h1, h2, J = fourier(R, b)
            t2 = b * b / r2
            Jn = J / (b ** 4 / r2 ** 3)
            hn = h1 / (b ** 2 / r2 ** 2)
            # 3. exact identities
            assert h1 == h2, "h1 != h2"
            assert fourier(R, -b)[3] == J and fourier(R, -b)[1] == h1, "F not even in b"
            if Rint[0] == 0 or Rint[1] == 0:
                assert h1 == 0, "axis field not exactly zero"
            # 2. leading coefficients (relative error must shrink like t^2)
            assert abs(Jn - (384 * s2 ** 2 - 160)) < 4000 * t2, "J leading term"
            assert abs(hn - 16 * s2) < 200 * t2, "h leading term"
            # 4. remainder in units of t^6 r^-2
            RJ = J - leading_J(R) * b ** 4
            Rh = h1 - leading_h(R, b)
            rem6 = (abs(RJ) + abs(Rh)) / (t2 ** 3 / r2)
            max_rem6 = max(max_rem6, rem6)
            print(f"  R={Rint} b={float(b):.0e}: J/(b^4 r^-6)={float(Jn):.6f}"
                  f" (target {float(384 * s2 ** 2 - 160):.6f}),"
                  f" h/(b^2 r^-4)={float(hn):.6f} (target {float(16 * s2):.6f}),"
                  f" remainder/(t^6 r^-2)={float(rem6):.1f}")
    print(f"  max remainder coefficient in units of t^6 r^-2: {float(max_rem6):.1f}  (< 4000)")
    assert max_rem6 < 4000
    # For t <= 1e-8 the remainder is therefore <= 4000 t^6 <= 4e-5 t^5 << 1e6 t^5.

    print("\n== 5. Internal field coefficient (23), (24), (32) ==")
    b, eta = Fr(1, 1000), Fr(1, 10 ** 9)
    B = (1 / (8 * (b + eta) ** 2) - 1 / (8 * b ** 2)) / 2
    beta = -Fr(1, 8) * b ** -3
    quad = abs(B - beta * eta) / (b ** -4 * eta ** 2)
    print(f"  B(eta)={float(B):.9e}  beta*eta={float(beta * eta):.9e}"
          f"  |B - beta eta| / (b^-4 eta^2) = {float(quad):.4f}  (bound 1)")
    assert quad < 1

    print("\n== 6. Sensitivity to eta, eq. (33) (bound 100) ==")
    b, eta = Fr(1, 10000), Fr(1, 10 ** 12)
    R = (Fr(1), Fr(1))
    c0, c1 = fourier(R, b), fourier(R, b, eta, Fr(0))
    sens = [abs((c1[i] - c0[i]) / eta) for i in range(4)]
    print("  d(A,h1,h2,J)/d eta =", [f"{float(s):.4f}" for s in sens])
    assert max(sens) < 100

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
