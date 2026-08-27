# Independent soundness audit of arXiv:2608.23506v2

**Auditor:** Claude Fable 5 (Anthropic), operating in Claude Cowork on 27 August 2026, at the request of the author.
**Scope:** internal logical soundness of the proof of Theorem 1.1 and independent recomputation of every constant that can be recomputed; not a peer review, and not a verification of the cited external results (NP-completeness of Independent Set on planar cubic graphs; existence of polynomial-area orthogonal grid drawings).
**Reproducibility:** every quantitative claim below is checked by `riesz2d_s2_exact_multipole_check.py` and `riesz2d_s2_exact_lattice_flip_check.py` in this directory (Python standard library only; exact rational arithmetic). Their outputs are stored alongside as `.txt` files.
**Note on provenance:** this document was written by an AI system. The author reviewed it and decided which recommendations to adopt in version 3 of the paper; see `../paper/CHANGES_v2_to_v3.md` and the *Status after revision* section at the end.

---


*Riesz Energy Subset Selection in the Euclidean Plane Is NP-Hard: A Reduction from the Ising Model on Planar Cubic Graphs* (M. Emmerich, v2, 25 Aug 2026)

Reviewed against the PDF, the TeX source, and the companion repository. Every quantitative claim that could be recomputed independently was recomputed in exact rational arithmetic (Python `fractions`), which at s = 2 is a strictly stronger check than the high-precision floats used in the repository.

## Verdict

I did not find a flaw that breaks the proof. The architecture (selector forcing → tree normalization by component flips → exact Fourier bookkeeping → threshold with a fixed source gap) is logically closed, the central algebraic identity is correct, and every numerical inequality in the hierarchy holds, usually with several orders of magnitude of slack. What I found is one arithmetic slip in a remark (harmless), one place where the written justification is looser than the numbers it certifies (the objective-box U-turn in the component-flip bound), and a handful of presentational points. None of these changes the truth of Theorem 1.1 as stated. Details below, ordered by importance.

## What was verified

**Source reduction (Prop. 3.1).** With σ = 2x − 1 and G cubic, Σ_E σ_uσ_v = 4q − 6t + 3n/2 and Σ_V σ_v = 2t − n, giving H_B = n/2 + 4(q − t) exactly as in (2). The gap of 4 between α ≥ K and α ≤ K − 1 in (4)–(5) follows. H_0 is an integer since n is even.

**Selector energetics (§4.1).** I_2 = b⁻²/8, side pair = b⁻²/4, I_3 = 5b⁻²/8, I_4 = 5b⁻²/4; the four transfer gains 3/8, 1/2, 1/2, 5/8 are correct and (3,1)→(2,2) is indeed the minimum. The argument that overfull cells with 3 or 4 points have internal energy exactly I_j (all 3-corner choices are congruent; 4 is unique) and that the best removal/insertion yields exactly I_{j−1}, I_{j'+1} is right, so the transfer gains are exact rather than lower bounds.

**Forcing lemma (Lemma 4.1).** The adverse cross-cell bound (2M − 1)(1 − 4b)⁻² < 3M is valid because every point lies within √2·b(1 + 10⁻⁶) of its centre, and the perturbation changes the six internal pair energies by roughly 3·10⁻⁶ b⁻², far below the 10⁻³ b⁻² allowance. Termination and the resulting normal form are sound. With b = (10⁸M³)⁻¹ the forcing condition (13) holds with ratio 10¹⁶M⁵/64.

**Multipole law (Lemma 5.1, Appendix B).** This is the load-bearing algebra, so I checked it two ways. By hand: the Fourier projections of c², −12a²c, 16a⁴ give 32, −192(x²+y²), 1536x²y² respectively, so [t⁴]J/r⁻² = 1536x²y² − 160 = 384 sin²2θ − 160; the t² term of J vanishes because Σ_{ε,η} εη = 0 and Σ_σ σ = 0; the field terms give 32xy and 256xy. Exactly: evaluating (14)–(16) in rational arithmetic at b = 10⁻³ and 10⁻⁴ for R ∈ {(1,0),(0,1),(1,1),(2,1),(3,1),(1,−1),(5,2),(7,3)} reproduces −160, +224, 85.76, −21.76, … to the expected order, and h/(b²r⁻⁴) = 16 sin 2θ. Three structural facts the proof relies on hold *exactly*, not just to leading order: F_R(σ,τ) = F_R(τ,σ) so h⁽¹⁾ = h⁽²⁾; F_R is even in b so all odd orders vanish; for axis-aligned R the field is identically zero. The Figure 3 values (J = −0.035776, +0.006145; F = 4.201188, 4.272740, 2.190833, 2.055895 at b = 0.12) are reproduced to all printed digits.

**Remainder constant (19).** Because F_R is even in b, the true remainder after the t⁴ term is O(t⁶ r⁻²). The measured coefficient is ≤ 3712 (attained at 45°), so |R_J| + |R_h| ≤ 4·10³ t⁶ r⁻² ≤ 4·10⁻⁵ t⁵ r⁻² at t ≤ 10⁻⁸. The stated 10⁶ t⁵ r⁻⁷·b⁵ bound is therefore true with about ten orders of magnitude to spare, and (20)–(22) follow.

**Field compensation (§6, Lemma 6.1).** B(η) = (1/16)[(b+η)⁻² − b⁻²] and β = −b⁻³/8 are correct; the measured constant in (32) is ≈ 0.19 b⁻⁴η² (bound: 1). The measured sensitivity of a two-spin Fourier coefficient to η is ≈ 10⁻³ per unit η at R = (1,1), b = 10⁻⁴ (bound in (33): 100). The decomposition G̃_i − g_i = [B(η_i) − βη_i] + Σ_j Δh_ij is exact, and both totals are below λ/100 by factors of order 10⁶ or more under (29). The identity η* = (g − G)/β is rational; bit lengths are polynomial.

**Objective box (Lemma 8.1).** Σ_{k≥2} k·k⁻⁶ = ζ(5) − 1 = 0.03693, so the halo bound 8.31 b⁴ and 19 < λ/b⁴ < 37 hold. The actual value is λ ≈ 29.18 b⁴ (D → ∞), and the exact rational computation at D = 2, b = 1/500 gives 29.270918 b⁴, matching the table in Appendix A.

**Component flips (Lemma 9.1).** I recomputed the absolute cross-cut coupling load (leading-order J, and separately the crude 225 r⁻⁶ bound the paper uses) for every cut type: straight, root junction with two perpendicular branches (crude 70.2 b⁴), root junction with one collinear branch (41.5 b⁴), and all cuts inside the objective-box template. The worst same-tree load is the cut on the unit step at the end of the q-prefix: 62.9 b⁴ crude (58.3 actual), plus 32.7 b⁴ cross-tree. Every cut stays below the paper's 95 + 37 budget, and the favorable edge is 160 b⁴, so the margin claim > 25 b⁴ is correct.

**Remote terms (46), normal form (Prop. 10.1), threshold (§11), NP membership (§12).** Arithmetic in (46) checks (2.25·10⁻⁴ L̄⁻⁶ b⁴). The classification (C1)–(C5) is exhaustive given the exact Fourier decomposition, and E_0 (perturbed constants plus perturbed same-tree couplings) is an exactly computable rational. The two threshold inequalities have slack 1.95λ and 1.95λ respectively. The verifier argument is fine.

## Issues found

### 1. Arithmetic slip in the "elementary certification" remark (p. 19) — harmless

The remark states Σ_{k=2}^{5} k⁻⁵ + ∫_5^∞ x⁻⁵ dx < 0.037. The left side is 0.036642 + 0.000400 = 0.037062, which is *not* below 0.037. The underlying constant ζ(5) − 1 = 0.036928 is fine, and even with 0.03706 the budget (44) evaluates to 94.28 b⁴ < 95 b⁴, so nothing downstream changes. Fix: extend the finite sum to k = 6 and integrate from 6 (gives 0.036983 < 0.037), or simply cite ζ(5) − 1.

### 2. The objective-box U-turn is not covered by "the ordinary bend estimate" as written — numbers hold, justification should be made explicit

The remark after Figure 4 says the reconnecting unit step "inside its own consistency tree … is covered by the ordinary bend estimate," and §9.1 defines local pairs as belonging to a "straight/bend/junction neighborhood" or an objective-box prefix pair. The q-prefix (on y = 1) and the continuation ray (on y = 0, x ≥ 1 + D) form two parallel arms at distance 1 joined by a unit segment: two bends, not one, and the two arms are non-adjacent segments at distance 1, not ≥ D. So these pairs are neither "remote" in the sense of (46) nor a single right-angle sector. In particular the pairs (q_{D−1}, (1+D, 0)) and (q_D, (2+D, 0)) are *same-tree 45° pairs* with antiferromagnetic coupling +28 b⁴ each — larger than the entire single-sector budget 225·Q_6 ≈ 41 b⁴ that "ordinary bend" would suggest.

The relevant lattice sum is Σ_{k≥1}(k+1)(k²+1)⁻³ ≈ 0.2796, i.e. ≈ 62.9 b⁴ with the 225 constant, and since every other feature is ≥ D away from this region, a cut anywhere in the U-turn sees at most 0.037 + 0.2796 = 0.317 < 0.419, so (44)/(45) remain valid. I recommend adding this sum explicitly to §9 (or to the Figure 4 remark) and adding the U-turn pairs to the definition of "local" in §9.1, so that the reader does not have to discover that the extra sectors in (44) happen to absorb it.

Relatedly, the sentence "a degree-three root can expose at most two perpendicular sectors and one opposite-collinear sector in addition to the straight tail" over-counts: three of four axis directions give, relative to the cut branch, either two perpendicular branches or one perpendicular plus one collinear, never all three. This is conservative and harmless, but it is exactly the slack that currently covers item 2, so it is worth saying so.

### 3. Lemma 6.1(iii): a factor 2 is silently absorbed

(33) bounds the change of *any single* Fourier coefficient by 100(|η_i| + |η_j|), but each pair contributes to two field coefficients (one per endpoint), so Σ_i |Σ_j Δh_ij| ≤ 2·100·Σ_{i<j}(|η_i|+|η_j|) = 30000 M³b⁵. Still smaller than λ/100 by a factor ~ 6·10⁴/M³·… under (29), so no consequence, but the displayed line "100 Σ_{i<j}(…) < 15000 M³b⁵" should either carry the 2 or say that 100 already absorbs both coefficients (the measured constant is ≈ 12, so either reading is safe).

### 4. Remainder bound (19) is asserted rather than derived

The proof of Lemma 5.1 says the tail is "bounded uniformly by the round constant 10⁶ t⁵ after the sixteen point-pair terms and the Fourier factor 1/4 are included" without showing the computation. Since it is true by a factor of ~10¹⁰, the simplest repair is to observe that F_R is even in b (substitute (ε,η) → (−ε,−η)), so the remainder is genuinely O(t⁶ r⁻²), and then give the one-line geometric-series bound Σ_{n≥5} (6t)ⁿ·16/4 ≤ 4·(6t)⁵/(1 − 6t) plus the finitely many degree ≥ 5 pieces of z², z³, z⁴ with |a| ≤ 2√2, c ≤ 8. That is a page of routine algebra and would make the lemma self-contained.

### 5. Presentation of the companion material

The two Lean files consist of `norm_num`/`ring`/`linarith` checks of the displayed arithmetic (e.g. 95 + 37 + 0.001 + 0.37 + 0.37 + 25 < 158.6) and a `native_decide` enumeration of Q₃; they do not formalize any lemma of the reduction. The Python script is a genuine finite-instance audit of the same architecture (it recomputes all fields and couplings from the rational coordinates and checks every tree edge), which is valuable. I would describe the Lean material in Appendix A as "arithmetic sanity checks" rather than leaving the reader to infer a formal verification. The symbolic-check script output confirms exactly the coefficients I verified independently.

### 6. Minor

- p. 10, Lemma 4.1: the hypothesis "cell centres mutually at distance ≥ 1" lives in the section preamble, not in the lemma statement; worth moving into the statement.
- §3: for NP-completeness of Independent Set on planar cubic graphs, a journal reference (Mohar 2001 is already cited, or Garey–Johnson–Stockmeyer 1976 for degree ≤ 3) is preferable to the Uehara technical report as the primary citation.
- Appendix A table: "maximum normalized remainder |R|/λ = 0.168" exceeds the theorem's 1/20; the text explains this is because b = 1/500 is far outside the asymptotic regime, but a reader skimming the table may take it as a counterexample. A one-line note in the table caption would help.
- The abstract's phrase "NP-complete" and the title's "NP-hard" are both correct; the Discussion already comments on this.

## Points examined and found sound (no action needed)

- Exactness of the Fourier decomposition (15)–(16) for the perturbed instance, hence exactness of E_0 and τ; only *bounds* use the multipole expansion.
- Invariance of λ under the 90° rotations and reflections of the template, since J depends on R only through sin²2θ and r while selector orientation is global.
- The leaf-side component of any tree-edge cut meets at most one objective box, so at most 37 b⁴ of cross-tree load is ever exposed.
- The root is never flipped, so the target field λ never enters the adverse budget; residual fields are bounded globally.
- Distinct non-adjacent segments of the scaled drawing are ≥ S apart (coordinates are multiples of S; planar, simple routing), so the "remote" classification is correct outside the objective-box template (see item 2 for inside it).
- Polynomial size: M ≤ 10⁶ L̄³ with L̄ = poly(n); all coordinates, b, η_i, E_0, τ have polynomial bit length.
- Direction of all sign conventions: negative axis coupling favours equal spins; positive 45° coupling favours opposite spins; +λσ_v for the root matches +Σσ_v in (1).

## Status after revision (v3)

All six items above were addressed by the author in version 3:

1. Certification remark corrected (sum through k = 6, integral from 6).
2. U-turn sum U₆ < 0.281 introduced as eq. (46); offset-sector bound (44) added; root case analysis corrected; "local pair" definition extended; Section 8 remark rewritten.
3. Factor 2 made explicit in Lemma 6.1.
4. Remainder bound of Lemma 5.1 derived explicitly (< 1.2·10⁵ t⁵ r⁻²), evenness in t noted.
5. Appendix A describes the Lean files as arithmetic sanity checks; remainder 0.168 λ of the finite example explained.
6. Lemma 4.1 hypothesis moved into the statement; Garey–Johnson–Stockmeyer citation added; finite-example caption note added.

The revised text was re-read by the auditor after the edits; the new constants (0.03698, 0.2794 + 0.0004, 94.5, 71.8, 7777, 6900, 1.2·10⁵, 30000 M³b⁵) are reproduced by the scripts in this directory and, as closed rational inequalities, by `../riesz2d_s2_audit_v3_constants.lean`.
