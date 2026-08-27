# Changes from arXiv v2 to v3

All changes follow the independent soundness audit (see `../audit/`). None of them changes a statement of a lemma or the theorem; they tighten the exposition and correct arithmetic in remarks. The full unified diff is in `tex_v2_to_v3.diff`.

## Corrections

- **Section 9, "elementary certification" remark.** The v2 bound `sum_{k=2}^{5} k^-5 + int_5^inf x^-5 dx < 0.037` is false (the left side is 0.03706). Replaced by the sum through k = 6 with the integral from 6, which gives 0.03698 < 0.037. The remark now says explicitly why the k = 5 version does not certify the constant. No downstream effect: ζ(5) − 1 = 0.03693 is what the proof uses.

## Gaps in justification made explicit

- **Section 9, objective-box U-turn.** The remark in Section 8 claimed that the unit step at the end of the q-prefix is "covered by the ordinary bend estimate" inside its own tree. It is not a single bend: the q-prefix and the corridor continuation are two parallel arms at distance 1, and the pairs (q_{D−1}, (1+D,0)) and (q_D, (2+D,0)) are same-tree 45° pairs with +28 b⁴ each. A new paragraph introduces the lattice sum U₆ = Σ_{k≥1}(k+1)(k²+1)^{-3} < 0.281 (new eq. (46)), shows that any cut in the U-turn region has local same-tree load < 72 b⁴, and the definition of "local pair" in Section 9.1 now lists the U-turn arms.
- **Section 9, position of the cut.** A new paragraph gives the offset sector bound Σ_{i≥a, j≥1}(i²+j²)^{-3} ≤ a^{-4} (new eq. (44)) and uses the ≥ 4D length of every straight run between features to show that at most one end of the cut run contributes; the root sentence ("two perpendicular sectors and one opposite-collinear sector") is replaced by the correct case analysis. The budget (45) gains a +0.001 term and still evaluates to 94.5 b⁴ < 95 b⁴.
- **Lemma 5.1, remainder bound (19).** The sentence asserting the constant 10⁶ t⁵ is replaced by an explicit derivation: geometric tail < 7777 t⁵, degree-5-to-8 pieces of z³ and z⁴ < 6900 t⁵, sixteen point-pair terms and the Fourier factor give |R_J| + |R_h| < 1.2·10⁵ t⁵ r^{-2}. The evenness of F_R in t (hence O(t⁶) true remainder) is noted. Appendix B updated to match.
- **Lemma 6.1 (iii).** The field-variation bound now carries the factor 2 for the two one-spin coefficients per pair (30000 M³b⁵ instead of 15000 M³b⁵), with the explicit inequality 22500 M³b⁶ + 30000 M³b⁵ < 4·10^{-4} b⁴ < λ/100 under (29).

## Presentation

- **Lemma 4.1** now states the hypothesis that selector centres are pairwise at distance ≥ 1 and k = 2M in the lemma itself.
- **Section 3** cites Garey–Johnson–Stockmeyer (1976) for planar max-degree-3 vertex cover alongside Mohar (2001) and Uehara (1996); new bibliography entry.
- **Section 8 remark** rewritten to say what the U-turn is and where it is charged.
- **Lemma 9.1 proof** references (47) and lists the cut types it covers.
- **Appendix A** describes the repository contents accurately (the Lean files are arithmetic sanity checks, not a formalization) and adds a paragraph explaining why the finite example's remainder 0.168 λ exceeds λ/20 without contradicting Proposition 10.1.
- **New Appendix C**: table of frequently used symbols with the section or equation where each is introduced.
- **Tool disclosure** mentions the Claude Fable 5 audit and points to the repository.
