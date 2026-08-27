# Support material and proof audit for arXiv:2608.23506

**Riesz Energy Subset Selection in the Euclidean Plane Is NP-Hard: A Reduction from the Ising Model on Planar Cubic Graphs**
Michael Emmerich, Faculty of Information Technology, University of Jyväskylä, Finland — August 2026

This repository accompanies the paper. It contains a finite worked instance of the reduction, symbolic and exact-rational cross-checks of the constants used in the proof, arithmetic sanity checks in Lean 4, and an independent AI soundness audit of version 2 together with the revised TeX source (version 3) that incorporates its recommendations.

## Layout

| Path | What it is |
|---|---|
| `cube8_riesz_reduction.py` | Builds the eight-vertex example of Appendix A (cube graph Q₃, 392 selector cells, b = 1/500, D = 2): orthogonal layout, selector trees, objective boxes, first-order field compensation with finite-decimal (hence rational) coordinates, exact two-state Hamiltonian at high precision, selector-forcing check, worst-case component-flip audit of every tree edge, and threshold separation over all 256 normalized states. Requires `mpmath`, `numpy` (optional `matplotlib`). Writes a JSON report, CSV tables, and a layout figure to `cube8_output/`. |
| `riesz2d_s2_symbolic_check.py` / `.txt` | SymPy expansion of the two-selector interaction: reproduces the exact coefficients 1536x²y² − 160 = 384 sin²2θ − 160, 32xy = 16 sin 2θ, 256xy = 128 sin 2θ of Lemma 5.1 and the ζ-constants of Section 9. Requires `sympy`. |
| `audit/claude_fable5_soundness_audit_v2.md` | Independent soundness audit of arXiv v2 performed with Claude Fable 5 (Anthropic) on 27 Aug 2026, with the list of items corrected in v3. Written by an AI system; see the provenance note inside. |
| `audit/riesz2d_s2_exact_multipole_check.py` / `.txt` | Exact-rational (`fractions.Fraction`) evaluation of the interaction (14) and its Fourier coefficients (16): Figure 3 values, leading coefficients of Lemma 5.1, the exact identities h⁽¹⁾ = h⁽²⁾ / evenness in b / zero axis field, the true O(t⁶) remainder (coefficient < 4000, versus the 10⁶ t⁵ budget of eq. (19)), and the constants of eqs. (23), (32), (33). Standard library only; asserts every claim. |
| `audit/riesz2d_s2_exact_lattice_flip_check.py` / `.txt` | ζ-constants, the corrected "elementary certification" arithmetic, the new U-turn sum U₆ < 0.281, λ ∈ (19, 37) b⁴ for several prefix lengths and the exact Appendix-A value 29.270918 b⁴, and the absolute cross-cut coupling load of Lemma 9.1 for every cut type (straight run, both root configurations, every cut inside the objective-box template including the unit-step U-turn). Standard library only; imports the multipole script; asserts every claim. |
| `riesz2d_s2_audit.lean` | Lean 4 / Mathlib arithmetic sanity checks of displayed inequalities of v2 (`norm_num`, `ring`, `linarith`). |
| `riesz2d_s2_audit_v3_constants.lean` | Same style, for the constants introduced or corrected in v3 (certification sums, U₆, budgets (45)/(47), remainder constants of Lemma 5.1, Lemma 6.1 slack). Closed rational inequalities decided by `norm_num`/`nlinarith`; written without a Lean toolchain at hand and not yet compiled — see below. |
| `riesz2d_s2_cube8_example.lean` | Exhaustive `native_decide` enumeration of Q₃: α = 4, min H_B = −12, next level −8, plus the finite-instance inequalities. |
| `paper/Emmerich_riesz2d_s2_planar_v3.tex` | Revised TeX source (v3). |
| `paper/tex_v2_to_v3.diff` | Unified diff against the arXiv v2 source. |
| `paper/CHANGES_v2_to_v3.md` | Human-readable change log. |
| `view.pdf` | Compiled paper (v2 at the time of writing). |

## What is and is not verified here

- The Python scripts in `audit/` verify, in exact arithmetic, every numerical constant of the general proof that can be evaluated independently of the source graph: the multipole coefficients and remainder, the lattice sums, the objective coefficient, and the component-flip budgets for all cut geometries. They do not verify the combinatorial routing lemma (Lemma 7.1) or the cited NP-completeness of Independent Set on planar cubic graphs.
- `cube8_riesz_reduction.py` verifies the complete architecture on one finite instance with human-scale constants (b = 1/500, D = 2). Its remainder 0.168 λ exceeds the asymptotic bound λ/20 of Proposition 10.1 because the instance is far outside the asymptotic regime; the relevant check there is the explicit threshold separation, which holds.
- The Lean files are **arithmetic sanity checks**, not a formalization of the reduction. Each theorem is a closed rational inequality or a finite enumeration.

## Running

```bash
# exact-arithmetic audit (no dependencies)
cd audit
python3 riesz2d_s2_exact_multipole_check.py
python3 riesz2d_s2_exact_lattice_flip_check.py

# finite instance (needs mpmath, numpy; matplotlib optional)
python3 cube8_riesz_reduction.py --outdir cube8_output

# symbolic coefficients (needs sympy)
python3 riesz2d_s2_symbolic_check.py

# Lean (needs a Mathlib project; copy the files into it)
lake env lean riesz2d_s2_audit.lean
lake env lean riesz2d_s2_audit_v3_constants.lean
lake env lean riesz2d_s2_cube8_example.lean
```

## Attribution

The audit in `audit/` and the v3 revisions to the TeX source were produced with Claude Fable 5 (Anthropic) in Claude Cowork on 27 August 2026 at the request of the author, who reviewed them and retains full responsibility for the paper. Earlier development used ChatGPT (OpenAI) for adversarial proof audits and symbolic cross-checks, as stated in the paper's tool disclosure.
