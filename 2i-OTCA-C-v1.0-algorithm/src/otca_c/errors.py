"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: none (defines exception types only; see the
module that raises a given exception for the specific result -- a
numbered theorem of the manuscript, or a supplementary repository
result -- whose domain is being enforced).

Exception types for 2i-OTCA-C. Kept in their own module so that callers can
catch them specifically without importing the rest of the package.
"""


class DomainError(ValueError):
    """Raised when the requested inputs fall outside the domain on which
    this component's underlying theorems (existence, uniqueness, global
    monotonicity -- see `circular/existence_uniqueness_minimizer/` and
    `circular/global_monotonicity/` at the repository root) are proven.

    This is raised instead of silently extrapolating, in keeping with this
    project's standing rule that a result is only ever reported within the
    region it has actually been established for.
    """


class ConvergenceError(RuntimeError):
    """Raised when the numerical root-finder fails to locate the minimizing
    semi-latus rectum `p*` within the safety-checked bracket, despite the
    input falling inside the proven domain of existence and uniqueness
    (the repository's supplementary existence/uniqueness result -- see
    `circular/existence_uniqueness_minimizer/` and `STATUS.md`; not a numbered
    lemma of the manuscript). This should not occur for well-posed inputs; if it does, it
    indicates a numerical edge case (e.g. extreme radius ratio) worth
    reporting rather than a silent wrong answer.
    """
