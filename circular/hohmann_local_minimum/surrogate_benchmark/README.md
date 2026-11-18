# Benchmark: quadratic curvature surrogate vs. exact cost

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

## What this folder contains

`benchmark_quadratic_surrogate.py` is a standalone numerical benchmark
supporting the discussion of practical relevance in Section 7.2 of the
published article. It quantifies,
rather than merely asserts, the claim that Theorem 2's closed-form
curvature coefficient $K(r_1,r_2)$ gives a cheap and accurate surrogate
for the exact two-impulse transfer cost over a bounded range of angular
deviations from Hohmann.

## The test problem

By Theorem 1, $V_R(\varphi)$ is strictly increasing on $(0,\pi)$, so the
*unconstrained* problem $\min_\varphi V_R(\varphi)$ always has its minimum
at $\varphi=0$ (Hohmann) — for both the exact cost and the surrogate.
Comparing the two methods on that problem would be trivial and
uninformative. The script instead uses the *constrained* problem

$$
\min_\varphi V_R(\varphi) \quad\text{subject to}\quad \varphi \ge \varphi_{\min},
$$

for a representative ratio $R=1.5$ and $\varphi_{\min}\in\{5^\circ,10^\circ,
15^\circ,20^\circ,25^\circ,30^\circ\}$ — modeling, for instance, a
phasing or rendezvous-window constraint that prevents the transfer angle
from reaching the Hohmann value exactly. Since both the exact cost and the
surrogate are confirmed strictly increasing on this range (checked
directly in the script, not assumed), the constrained optimum of both
coincides with $\varphi=\varphi_{\min}$, which makes the comparison
well-defined and non-trivial.

## What is measured

- **Accuracy**: the relative error of the surrogate $V_R(0)+K\varphi^2$
  against the exact cost (obtained by numerically minimizing the full,
  unexpanded two-impulse cost over the semi-latus rectum $p$), evaluated
  at $\varphi=\varphi_{\min}$ for each tested $\varphi_{\min}$.
- **Speed**: a timing comparison (2000 repeated evaluations at
  $\varphi=15^\circ$) between the exact (numerically minimized) cost and
  the surrogate (a single algebraic evaluation).
- **Agreement of the constrained optima**: confirmed by directly checking
  monotonicity of both cost functions on each tested range, rather than
  assumed from Theorem 1 alone (the surrogate's monotonicity in $\varphi>0$
  is immediate, but is checked numerically here alongside the exact cost
  for a like-for-like comparison).

## Status

- **Scientific status**: a numerical benchmark, not a proof; it supports
  (does not replace) the closed-form guarantee of Theorem 2.
- **Reproducibility status**: `REPRODUCED`. Running the script gives a
  surrogate relative error from $0.0003\%$ ($\varphi_{\min}=5^\circ$) to
  $0.24\%$ ($\varphi_{\min}=30^\circ$) at $R=1.5$, growing as predicted by
  the $O(\varphi^4)$ remainder of Theorem 2's expansion, and a
  per-evaluation speedup on the order of several hundred times (observed
  range $578\times$-$907\times$ across repeated runs and tested
  $\varphi_{\min}$ values on the reference machine used here) for the
  surrogate over the exact, numerically minimized cost. The relative-error
  values are deterministic (a numerical minimization, not a timing); the
  speedup fluctuates with system load and is reported only as an order of
  magnitude.

## Scope note

This is a benchmark of the surrogate's accuracy and speed on one
representative test problem (two-impulse, circular-to-circular,
$R=1.5$). It is not a claim that this surrogate has been deployed inside
any specific external mission-design tool, optimizer, or software
package; the manuscript's Section 7.2 states this scope explicitly. The
reported speedup factor is implementation- and hardware-dependent and
should be read as an order of magnitude, not a universal constant.

## Reproduction

Requires `numpy` and `scipy`. Run directly:

```
python3 benchmark_quadratic_surrogate.py
```

Runs in a few seconds.
