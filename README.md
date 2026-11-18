# Two-Impulse Transfer Cost Between Circular Orbits

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Interactive demonstrator:** [`2i-OTCA-C-v1.0-demo/index.html`](2i-OTCA-C-v1.0-demo/index.html)
-- open this file in a browser to compute transfers interactively, review
the theoretical results, and read the Algorithms documentation, without
installing anything.

This repository accompanies the paper *"Monotonicity, Convexity, and
Curvature of the Optimal Two-Impulse Transfer Cost Between Circular Orbits"*
(part of the 2i-OTC-C Study), archived separately (see `CITATION.cff`). It
contains the mathematical proofs, their corresponding computational
certificates, a reference software implementation
(`2i-OTCA-C-v1.0-algorithm/`), and an interactive demonstrator
(`2i-OTCA-C-v1.0-demo/`) for two-impulse transfers between circular orbits.

See `STATUS.md` for the scientific and reproducibility status of every
theorem, lemma, and conjecture in this repository, including which results
are the article's numbered theorems and which are supplementary results
that are not part of the published manuscript's numbering. The `Numbering
note` section of `STATUS.md` provides the full explanation.

## Structure

- `circular/` -- proofs and computational certificates for two-impulse
  transfers between circular orbits: existence and uniqueness of the
  minimizer, global monotonicity, universal convexity, and the
  curvature/inflection structure (including the still-open extension to
  general R).

  **Note:** Two of this folder's internal theorem labels ("Theorem 1" /
  "Theorem 2" in `hohmann_local_minimum/` and `global_monotonicity/`) are
  swapped relative to the article's numbering. In addition, two other
  folders (`existence_uniqueness_minimizer/` and
  `curvature_sign_change_general_R/`) contain results that are not numbered
  theorems or conjectures in the published article. See `STATUS.md` and
  `docs/mathematical_conventions.md` before citing a folder name as "the
  article's Theorem N."

- `2i-OTCA-C-v1.0-algorithm/` -- **2i-OTCA-C** (*Two-Impulse Orbital
  Transfer Cost Algorithm -- Circular Case*), a software component
  implementing the circular-case cost function, the admissible-angle
  reduction enabled by the proven monotonicity result, and a fast
  closed-form surrogate valid near the Hohmann transfer. Circular case
  only; see its own README for an explicit statement of scope and
  out-of-scope cases.

- `2i-OTCA-C-v1.0-demo/` -- interactive, bilingual (EN/FR) HTML
  demonstrator featuring live solutions, independent geometry/Lambert
  cross-checks, and a Keplerian-speed animation of the transfer arc.

- `docs/mathematical_conventions.md` -- the single source of truth for
  notation throughout the repository.

- `docs/proof_dependencies.md` -- the dependency map between theorems,
  lemmas, and the hypotheses on which they rely.

## Scope

This repository covers two-impulse, coplanar, impulsive-burn transfers
between circular orbits only. It does not model or optimize elliptical-orbit
transfers, orbital-plane changes, perturbations (e.g., J2), finite-duration
or finite-thrust propulsion, or launch-window, phasing, or rendezvous
constraints. It is a research software component, not a general-purpose
mission planner.

## License

This repository is dual-licensed:

- **Code and software** -- Research, Education and Non-Commercial Use
  License (NCSL-1.0). See `LICENSE-CODE.md`.
- **Theoretical and scientific content** (proofs, paper, Theory tab) --
  Creative Commons Attribution 4.0 International (CC BY 4.0). See
  `LICENSE-CONTENT.md`.

## Citation

See `CITATION.cff`.

> Speyser, F. G., & Vyatkina, K. (2026). *Monotonicity, Convexity, and
> Curvature of the Optimal Two-Impulse Transfer Cost Between Circular
> Orbits*. Zenodo. https://doi.org/10.5281/zenodo.24354529
