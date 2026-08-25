import Mathlib

namespace Riesz2DS2Audit

/-- The side pair exceeds the diagonal by exactly 1/8 in units of b^{-2}. -/
theorem selector_side_gap : (1 : ℚ) / 4 - 1 / 8 = 1 / 8 := by
  norm_num

/-- The worst (3,1) -> (2,2) balancing move gains 3/8 in units of b^{-2}. -/
theorem selector_balance_gap : (5 : ℚ) / 8 - 2 * (1 / 8) = 3 / 8 := by
  norm_num

/-- Exact reconstruction of the ++ two-spin value from Fourier coefficients. -/
theorem fourier_pp (pp pm mp mm : ℝ) :
    let A  := (pp + pm + mp + mm) / 4
    let h1 := (pp + pm - mp - mm) / 4
    let h2 := (pp - pm + mp - mm) / 4
    let J  := (pp - pm - mp + mm) / 4
    A + h1 + h2 + J = pp := by
  dsimp
  ring

/-- Exact reconstruction of the +- two-spin value. -/
theorem fourier_pm (pp pm mp mm : ℝ) :
    let A  := (pp + pm + mp + mm) / 4
    let h1 := (pp + pm - mp - mm) / 4
    let h2 := (pp - pm + mp - mm) / 4
    let J  := (pp - pm - mp + mm) / 4
    A + h1 - h2 - J = pm := by
  dsimp
  ring

/-- The conservative s=2 component-flip budget leaves more than 25 b^4. -/
theorem consistency_budget :
    (95 : ℚ) + 37 + 1/1000 + 37/100 + 37/100 + 25 < 793/5 := by
  norm_num

/-- Barahona source values are separated by four units. -/
theorem barahona_gap (H0 : ℝ) : H0 + 4 - H0 = 4 := by
  ring

/-- Exact-threshold completeness: a yes state plus residual fits below H0+2. -/
theorem threshold_completeness
    (A lam R : ℝ)
    (hlam : 0 < lam)
    (hR : R ≤ lam / 20) :
    A + R < A + 2 * lam := by
  linarith

/-- Exact-threshold soundness: a no state at H0+4 stays above H0+2. -/
theorem threshold_soundness
    (A lam R : ℝ)
    (hlam : 0 < lam)
    (hR : -lam / 20 ≤ R) :
    A + 2 * lam < A + 4 * lam + R := by
  linarith

end Riesz2DS2Audit
