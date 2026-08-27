import Mathlib

/-!
Arithmetic sanity checks for the constants introduced or corrected in
version 3 of arXiv:2608.23506 (Section 9 and Lemmas 5.1, 6.1).

These are decidable rational inequalities discharged by `norm_num` /
`linarith`; they certify the displayed arithmetic, not the reduction itself.
Added after the independent soundness audit with Claude Fable 5 (Anthropic),
27 August 2026.  The file was written without access to a Lean toolchain;
every statement is a closed rational inequality that `norm_num` decides, but
it has not been compiled by the author of the audit.
-/

namespace Riesz2DS2AuditV3

/-- The v2 "elementary certification" of the straight-tail constant fails:
    sum_{k=2}^{5} k^-5 + int_5^inf x^-5 dx = 0.03706 is NOT below 0.037. -/
theorem straight_tail_v2_certification_fails :
    (1/32 + 1/243 + 1/1024 + 1/3125 : ℚ) + 1/2500 > 37/1000 := by
  norm_num

/-- The corrected certification: sum_{k=2}^{6} k^-5 + int_6^inf x^-5 dx < 0.037. -/
theorem straight_tail_certification :
    (1/32 + 1/243 + 1/1024 + 1/3125 + 1/7776 : ℚ) + 1/5184 < 37/1000 := by
  norm_num

/-- Partial U-turn sum: sum_{k=1}^{5} (k+1)(k^2+1)^-3 < 0.2794. -/
theorem uturn_partial_sum :
    (2/8 + 3/125 + 4/1000 + 5/4913 + 6/17576 : ℚ) < 2794/10000 := by
  norm_num

/-- U-turn tail: sum_{k>=6} (k^-5 + k^-6) <= (6^-5 + int_6^inf x^-5) + (6^-6 + int_6^inf x^-6) < 0.0004. -/
theorem uturn_tail_bound :
    (1/7776 + 1/5184 : ℚ) + (1/46656 + 1/38880) < 4/10000 := by
  norm_num

/-- Hence U_6 < 0.281. -/
theorem uturn_constant : (2794/10000 : ℚ) + 4/10000 < 281/1000 := by
  norm_num

/-- Local same-tree budget (45): 225 (0.037 + 2*0.181 + 0.020 + 0.001) < 95. -/
theorem local_tree_budget :
    (225 : ℚ) * (37/1000 + 2 * (181/1000) + 20/1000 + 1/1000) < 95 := by
  norm_num

/-- U-turn cut budget: 225 (0.037 + 0.281 + 0.001) < 72 < 95. -/
theorem uturn_cut_budget :
    (225 : ℚ) * (37/1000 + 281/1000 + 1/1000) < 72 ∧ (72 : ℚ) < 95 := by
  constructor <;> norm_num

/-- The two same-tree 45-degree pairs of the U-turn: 2 * 224/8 = 56 b^4. -/
theorem uturn_diagonal_pairs : (2 : ℚ) * (224 / 8) = 56 := by
  norm_num

/-- Sector offset (44): (3 pi / 16)(a^-5 + a^-4/4) <= a^-4 for a >= 2, using pi < 3.15.
    Stated with the rational majorant p of pi and the worst case a = 2. -/
theorem sector_offset_constant (p : ℚ) (hp : p ≤ 315/100) (hp0 : 0 ≤ p) :
    3 * p / 16 * (1/2 + 1/4) ≤ 1 := by
  nlinarith

/-- Far end of a long run: (2D)^-4 < 10^-3 already for D >= 3. -/
theorem far_end_negligible : ((2 * 3 : ℚ))⁻¹ ^ 4 < 1/1000 := by
  norm_num

/-- Remainder bound of Lemma 5.1: geometric tail (6t)^5/(1-6t) < 7777 t^5 at t <= 1e-8,
    i.e. 6^5 / (1 - 6e-8) < 7777. -/
theorem geometric_tail_constant : (6 : ℚ) ^ 5 / (1 - 6 / 10 ^ 8) < 7777 := by
  norm_num

/-- Degree-five-to-eight pieces of z^3 and z^4 with |a| <= 2 sqrt 2, c <= 8, sqrt 2 < 1.415,
    t <= 1e-8: below 6900 t^5. -/
theorem high_degree_pieces (s : ℚ) (hs : s ≤ 1415/1000) (hs0 : 0 ≤ s) :
    6 * (2 * s) * 64 + 32 * (16 * s) * 8
      + (512 + 12288) / 10 ^ 8 + 8 * (2 * s) * 512 / 10 ^ 16 + 4096 / 10 ^ 24 < 6900 := by
  nlinarith

/-- Sixteen point-pair terms, Fourier factor 1/4, two coefficients:
    2 * 4 * (7777 + 6900) < 1.2e5 <= 1e6. -/
theorem remainder_total : (2 : ℚ) * 4 * (7777 + 6900) < 120000 ∧ (120000 : ℚ) ≤ 10 ^ 6 := by
  constructor <;> norm_num

/-- Lemma 6.1 with the factor 2 for the two field coefficients per pair:
    22500 M^3 b^6 + 30000 M^3 b^5 < 4e-4 b^4 < 0.19 b^4 when M^3 b <= 1e-8 and M^3 b^2 <= 1e-16. -/
theorem field_compensation_slack (x y b4 : ℚ) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (hx' : x ≤ 1 / 10 ^ 16) (hy' : y ≤ 1 / 10 ^ 8) (hb : 0 < b4) :
    22500 * x * b4 + 30000 * y * b4 < 19/100 * b4 := by
  nlinarith

end Riesz2DS2AuditV3
