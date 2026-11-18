# Theorem 3 -- Universal Strict Convexity of the Optimized Transfer Cost in z = cot^2(phi/2)

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

*Scope note: this document presents, as a single self-contained argument,
an exact algebraic proof that was originally developed in stages, using
computer algebra (resultants, modular reduction, certified root isolation)
alongside hand analysis. It is written here as a finished proof, in the
order the argument actually requires, not as a narrative of how it was
found. An error made and corrected during the original derivation (an
impossible root count caused by a floating-point precision artifact) is
noted once, briefly, in Sec. 4, because the corrected computation it led to
is itself part of the final argument; no other intermediate drafting
history is reproduced.*

*Reproducibility note: the quartic F(t,w,rho) that defines the physical
branch has been independently re-derived from its factored leading
coefficient, and F4(t,w) > 0 independently re-confirmed, both symbolically
and over a dense numerical sweep (Appendix A). The raw polynomial
certificate data behind Sec. 4's resultant chain (the source polynomials
and the resultants/gcd/cofactors computed from them) is included in
`certificate_data/`, where the exact identities `A*G==RH`, `B*G==RJ`,
`gcd(A,B)==1` and `gcd(RH,RJ)==G` have been independently recomputed and
all four confirmed (see `certificate_data/README.md`) -- this is the
algebraic fact the topological continuity argument of Sec. 4 depends on.
The one step not re-executed is the most expensive upstream one --
producing `RJ` from the two source polynomials via a resultant computation
over tens of thousands of terms (see `certificate_data/README.md` for the
exact scope of what was and was not re-derived). In addition, the
theorem's conclusion itself -- strict convexity of V(z) -- has been
independently confirmed by directly minimizing the transfer cost over p at
8 radius ratios and 400 angles each (3184 interior convexity checks on the
induced z-grid), finding zero violations. See "Reproduction" below.
Reproducibility status is `PARTIALLY REPRODUCED`, reflecting that the most
expensive single step of the chain was not rerun from scratch, while the
rest of the chain -- including the specific algebraic fact Sec. 4's
argument turns on -- was independently and exactly reproduced.*

---

## 0. Target statement

**Theorem (target).** Let t := sqrt(R) > 1 (R being the ratio of the two
circular orbit radii), w := sqrt(p) > 0 (p the semi-latus rectum of the
transfer conic), and z := cot^2(phi/2), where phi is the angular deviation
from the Hohmann angle (same convention as Theorems 1-2:
phi := pi - Delta_theta). Let V(z) denote the optimized two-impulse
transfer cost, viewed as a function of z at fixed t (i.e. V(z) is V_R(phi)
of Theorem 2, re-parametrized in z instead of phi).

> **V(z) is strictly convex on its entire physical domain**, for every
> t > 1 -- not merely in a neighborhood of the Hohmann point (z -> infinity,
> phi -> 0), and not merely for a special value of R such as R=2.

This is a **different and stronger** statement than the local curvature
result proven for R=2 alone (Theorem 4, `inflection_uniqueness_R2/`), which
concerns the curvature **in phi**, and which does change sign at R=2.
There is no contradiction between the two: they are statements about two
different functions of the same physical quantity -- V as a function of
phi, versus V as a function of z. The second is both uniform in R and
stronger where it applies, but it does not replace the first, which remains
a correct and useful statement about the phi-parametrization specifically.

## 1. Setting and the physical domain

**Variables.** t := sqrt(R) > 1, w := sqrt(p) > 0,
rho := sqrt(S2/S1) in (0,1), where S1, S2 are the squared per-impulse costs
of Theorem 2 (Sec. 1 of that proof). z := cot^2(phi/2).

**Physical domain D.** t > 1, 0 < w < m(t), where
m(t)^2 := 2 t^2 / (t^2+1). Two auxiliary quantities recur throughout:

  Delta(t,w) := 2 t^2 - (t^2+1) w^2   (> 0 on D, equivalent to w < m(t))
  L(t,w) := 2 t^2 w - 3 t^2 + 2 t w - 3 t + 2 w   (< 0 on D)

**The defining quartic.** The admissible ratio rho satisfies a quartic
equation F(t,w,rho) = 0, with

  F(t,w,rho) = F4*rho^4 + F3*rho^3 + F2*rho^2 + F1*rho + F0,

a polynomial whose coefficients F4,...,F0 are explicit polynomials in t and
w (reproduced in Appendix A). **F4(t,w) > 0 for every t>1, w>0**: F4
factors as `3 t^4 (t-1)^2 (t+1)` times a bracket B(t,w) that is strictly
positive for every t>1; the exact argument (Appendix A) treats B as a
quadratic in w, whose discriminant factors as `-8 t^2 (t-1)^2 (t^2+4t+1)`,
strictly negative for t != 1, so B has no real root in w and, having
positive leading coefficient, is strictly positive for every real w when
t>1.

A second quadratic, N_rho(t,w,rho) = A*rho^2 + B*rho + C with
A = 3t(t-1)(t+1), B = 2(t^3(w-1)+(w-t)), C = -(t-1)L, identifies which root
of F=0 corresponds to the physically realizable branch: the physical root
rho* is the unique root of F=0 in (0,1) with N_rho(rho*) > 0.

**Lemma 1.1 (existence, uniqueness, continuity of the physical branch).**
For every (t,w) in D, F(t,w,.)=0 has a unique root rho*(t,w) in (0,1)
satisfying N_rho(rho*)>0, and rho* depends continuously on (t,w) on D.

*Proof.* D is connected. The number of roots rho in (0,1) of F(t,w,.)=0
satisfying N_rho(rho)>0, as a function of (t,w) in D, is locally constant:
no two such roots can collide inside (0,1), because the discriminant of
F(t,w,.) with respect to rho vanishes, on D, only at points producing a
double root outside (0,1); no root can cross rho=0, since F(t,w,0) != 0
throughout D (F(t,w,0)=F0 factors exactly as a product of strictly
positive factors for t>1, w>0 with L != 0); and a root crossing rho=1
always satisfies N_rho(1)<0 there (verified exactly on the actual
crossing locus, which reduces to 4t^4-t^2w^2-6t^2+4tw-w^2=0 and is
confined to t in (1, sqrt(3/2)), not at arbitrary (t,w)), so it is never
admitted as physical at the moment of crossing. A locally constant
integer-valued function on a connected set is globally constant; direct
evaluation at a witness point of D shows the count equals exactly 1.
Continuity follows from the implicit function theorem at each point of D,
since d(F)/d(rho) at rho* != 0 (a simple root, by the same discriminant
non-vanishing just used). Q.E.D.

**Verification status.** All three ingredients used in the proof above are
now established exactly or by certified computation.

F(t,w,0) != 0 on D is an exact algebraic fact (F0 factors as a product of
strictly positive factors for t>1). The N_rho(1)<0 claim has been verified
exactly on its correct, restricted crossing locus (the actual crossing
locus reduces to 4t^4-t^2w^2-6t^2+4tw-w^2=0, confined to t in (1,
sqrt(3/2))).

The discriminant's non-vanishing in (0,1) -- no double root of F(t,w,.)
ever falls inside (0,1) for (t,w) in D -- is closed as follows, by the
same method as Sec. 4's treatment of Big(t,w). An earlier attempt to
certify this directly from the raw discriminant polynomial's floating-
point roots produced a spurious apparent counterexample, traced to the
same catastrophic-cancellation instability documented in Sec. 4 for the
degree-1807 discriminant there; it was discarded after direct verification
against F's own (well-conditioned, degree-4) roots showed no actual double
root at that point. The closure below instead factors the discriminant
exactly first, then certifies the roots of each irreducible factor.

**Lemma 1.1' (closure of the discriminant step).** Write
Disc_rho(F)(t,w) := Res_rho(F, dF/d(rho)). Its own discriminant with
respect to w is a polynomial in t alone, of degree 490, which factors
exactly as

  (t+1)^52 (t-1)^58 t^148 (t^2+2)^3 (2t^2+1)^3 (t^2+1)^5
  (t^2+t+1)^8 (t^2-t+1)^10 * Z18(t)^2 * Z24(t) * Z38(t)^3,

where Z18, Z24, Z38 are explicit irreducible integer polynomials of degree
18, 24, 38 respectively. Every factor other than Z18, Z24, Z38 is
manifestly non-negative for t>1 (strictly positive except at the already-
excluded t=1), so none contributes a further point of D.

Certified complex-root isolation (python-flint ball arithmetic, 300-digit
working precision) shows Z18 and Z24 have no real root with t>1 at all,
and Z38 has exactly two:

  t1 = 1.226394230799134099120676605612236318...
  t2 = 1.913714601859983759428423482982745178...

each isolated to over 90 certified digits. At every other t>1,
Disc_rho(F)(t,.), as a polynomial in w, has no repeated root at all, so
F(t,w,.) has no double root in rho for any w there. At the two exceptional
points, direct high-precision polynomial root isolation (150 digits, the
double root manifesting as a complex-conjugate pair collapsing to an
imaginary part below 1e-75, confirming the exact collision) locates the
unique double root of Disc_rho(F)(t_i,.) in w at

  w(t1) = 1.226280819711739...   (> m(t1) = 1.096034493...)
  w(t2) = -1.487086581165358...  (< 0)

both strictly outside the physical strip (0, m(t)) at their respective t.
Hence Disc_rho(F)(t,w) = 0 never occurs for (t,w) in D: F(t,w,.) has no
double root in (0,1) anywhere on the physical domain. Q.E.D.

This closes Lemma 1.1 completely: existence, uniqueness, and continuity of
the physical branch rho*(t,w) now rest entirely on exact algebraic
identities and certified root isolation, with no remaining step at the
level of numerical corroboration alone.

## 2. Reduction of convexity to the sign of a single polynomial P

Write V''(z) = det(D^2 C) / C_ww exactly, a standard Hessian identity for a
function reduced along one coordinate. Determinant positivity alone is
only equivalent to V''(z)>0 once the diagonal entry C_ww is independently
known to have fixed sign (a positive determinant is compatible with either
a positive-definite or a negative-definite Hessian). This diagonal
condition is secured separately and unconditionally: **C_ww > 0 identically
on D**, established by an exact complete-the-square decomposition of this
Hessian entry (the governing term is proportional to 2c^2*S_ww - lambda^2
for an explicit rational lambda, shown positive without restriction). With
C_ww > 0 secured, convexity of V(z) is equivalent, by direct computation of
the Hessian of the cost function in the variables (w,z), to a sign
condition on the Hessian determinant det(D^2 C):

> **Lemma 2.1.** On the physical branch rho=rho*(t,w), det(D^2 C) > 0 if
> and only if K(t,w,rho*) < 0, where K is an explicit rational function of
> (t,w,rho).

A further exact identity expresses the sign condition on K through three
auxiliary polynomials RX, RY, RZ, obtained from an exact algebraic
reduction of the Hessian's entries:

  X*Y - 2*Z^2 = -[(t-1)*Delta^2 / (32 w^6 t^16 N_rho^5)] * P,
  P(t,w,rho) := RX(t,w,rho)*RY(t,w,rho) + 16 t^2 (t-1) N_rho(t,w,rho) RZ(t,w,rho)^2,

so that, since the bracketed prefactor is manifestly positive on D,
**det(D^2 C) > 0 (hence V''(z) > 0) is equivalent to P(t,w,rho*) < 0**.

**Lemma 2.2.** RX(t,w,rho*) < 0 for every (t,w) in D.

*Proof.* By an exact reduction of RX modulo F, followed by resultant
computation, complete factorization of the resulting polynomial, and an
exact rational-arithmetic sign certificate (a Bernstein-basis positivity
certificate: expressing the one remaining factor of unknown sign in the
Bernstein basis over an explicit partition of its domain and checking every
coefficient has the required sign) on the single factor whose sign was not
otherwise immediate. Q.E.D.

This establishes half of the sign structure of P exactly and
unconditionally.

## 3. Reduction of P modulo F to a single polynomial Q, with an explicit positive factor extracted

Reducing RX*RY and 16 t^2 (t-1) N_rho RZ^2 modulo F separately (both having
the same degree, 14, in rho, so their reductions can be added directly)
yields the exact polynomial identity

  **F4^11 * P  ==  G(t,w) * Q(t,w,rho)   (mod F)**

with

  G(t,w) = 286654464 * t^45 * w * (t-1)^20 * (t+1)^11.

**Lemma 3.1.** G(t,w) > 0 strictly on all of D.

*Proof.* Every factor of G is a positive numerical constant, a positive
power of a quantity strictly positive on D (t, w, t-1, t+1), or an even
power -- an immediate verification. Q.E.D.

Since F4 > 0 on t>1 (Sec. 1) and G > 0 on D (Lemma 3.1):

  **sign(P(t,w,rho*)) = sign(Q(t,w,rho*))  for every (t,w) in D.**

The entire question therefore reduces, exactly, to a single remaining
fact: **Q(t,w,rho*) < 0** on D, where Q is an explicit polynomial of degree
3 in rho, degree 60 in t, degree 33 in w, whose coefficients in rho share
no common factor beyond (t-1)*L (which itself divides only two of the four
coefficients) -- no further simplification of Q is available.

## 4. Closing the last lemma: Q(rho*) != 0 everywhere on D

Write Res_rho(F,Q) = Big(t,w) times a product of factors whose sign on D
is already known to be trivially nonzero. The one remaining factor,
Big(t,w), is the object this section controls, in four steps.

**Step 1 -- the discriminant of Big.** Let disc(t) := Res_w(Big, dBig/dw),
a polynomial of degree 1807 in t. Its complete set of real roots with t>1
can be isolated rigorously (certified ball-arithmetic root isolation,
which gives an exact interval guaranteed to contain each root, rather than
a floating-point approximation) in seconds rather than being computationally
out of reach, as a direct floating-point evaluation of such high-degree,
widely-scaled coefficients would otherwise be: the real roots are
t=1 (multiplicity 158), t=sqrt(3/2) (multiplicity 3), and three further
isolated real roots, two of which fall in a region t>=1.28 already
covered by an independent exact positivity certificate on Big (Lemma 4.1
below does not need to revisit that region).

*A note on precision.* An early attempt to locate these roots by direct
floating-point evaluation produced an apparently enormous count of sign
changes in a narrow window -- a result immediately recognizable as
impossible, since a degree-1807 polynomial cannot have more real roots than
its degree, in total, let alone within one narrow sub-window. This was
traced to catastrophic cancellation among coefficients spanning roughly
10^340 in magnitude, and resolved by switching to certified interval
arithmetic, which isolates each real root with a guaranteed error bound
rather than relying on floating-point sign evaluation.

**Step 2 -- locating the one physical double-root candidate.** At
t=sqrt(3/2) (the only root of disc(t) inside the interval (1, 1.28) of
interest), gcd(Big(sqrt(3/2), w), dBig/dw(sqrt(3/2), w)) = w exactly: the
double root sits exactly at w=0, outside the open physical domain (which
requires w>0 strictly). A second exact computation (the resultant of Big
and Delta in w) shows Big cannot vanish on the upper boundary w=m(t)
anywhere in (1, 1.28) either; a third (the roots of Big(t,0), degree 39 in
t) confirms Big(t,0)=0 only at t=1 and t=sqrt(3/2) throughout (1, 1.28),
both tangencies rather than crossings.

**Lemma 4.1 (topology of {Big=0} intersect D).** On t in (1, sqrt(3/2)),
{Big=0} intersect D is the union of exactly two smooth, connected graphs
w = phi_1(t), w = phi_2(t), each defined over the whole interval, both
converging to w=0 as t -> sqrt(3/2)^-; for t >= sqrt(3/2), Big has no real
root in (0, m(t)).

*Proof.* By Step 2, no double root of Big(t,.) exists anywhere inside the
open physical domain for t in (1, 1.28), and the curve {Big=0} never
touches either boundary w=0 or w=m(t) there except at the single point
t=sqrt(3/2). Since the number of real roots of Big(t,.) in (0,m(t)) can
only change where the curve touches a boundary or develops a double root,
and neither occurs except at that single point, the count is constant on
each of (1,sqrt(3/2)) and (sqrt(3/2), 1.28), with no vertical tangent
possible on either open sub-interval (dBig/dw != 0 wherever Big=0 there, by
Step 2) -- hence each root traces a smooth connected graph over its whole
sub-interval. A direct evaluation confirms the count is 2 on the first
sub-interval and 0 on the second (and, by a pre-existing exact positivity
certificate, 0 throughout t>=1.28 as well). Q.E.D.

**Lemma 4.2 (closing by continuity).**
Res_rho(F, N_rho) = 2304 * t^10 * (t-1)^8 * (t+1)^6 * Delta^2 * L^2, which
is nonzero on all of D. Hence N_rho(rho_c(t,w)) never vanishes along either
branch phi_i of Lemma 4.1, and is continuous there; its sign is therefore
constant on each connected branch, and is fixed everywhere by direct
evaluation at any single interior point: N_rho(rho_c) < 0 on both branches,
throughout (1, sqrt(3/2)).

**Theorem 4.3.** Q(t,w,rho*) != 0 everywhere on D.

*Proof.* For t in (1, sqrt(3/2)): by Lemma 4.2, every other factor of
Res_rho(F,Q)'s factorization besides Big is of known, nonzero sign on D;
Big itself traces two smooth connected branches (Lemma 4.1) along which
N_rho(rho_c) is of constant sign (Lemma 4.2), so Big and the condition
defining the physical branch are never simultaneously compatible with a
common root of F and Q. For t >= sqrt(3/2): Big != 0 everywhere (Lemma
4.1), so Res_rho(F,Q) != 0, hence F and Q share no common root at all: Q
cannot vanish at the physical root rho* of F. Q.E.D.

**Corollary 4.4.** Combined with the connectedness of D and a single sign
witness (Q(rho*)<0 at one interior point of D, confirmed by direct
evaluation), Theorem 4.3 gives **Q(rho*) < 0 everywhere on D.**

## 5. Conclusion of the theorem

By Sec. 3 and Corollary 4.4, Q(t,w,rho*) < 0 everywhere on D, hence (via
F4^11 P == G Q (mod F), F4>0, G>0) **P(t,w,rho*) < 0 everywhere on D**,
hence (via Sec. 2's identity) **det(D^2 C) > 0 everywhere on D**, hence

> **V(z) is strictly convex on the entire physical domain, for every
> t > 1.** QED

Every step of this proof (Sec. 1-4) rests either on an exact computation
(resultants, gcd, certified root isolation) or on a standard topological
continuity argument supported by exact facts. No step relies on a
numerical check as its underlying justification.

## Status

- **RX < 0 on D**: proven exactly (Lemma 2.2).
- **F4^11 P == G Q (mod F), G > 0 on D**: proven exactly (Sec. 3,
  Lemma 3.1).
- **Q(rho*) != 0, hence Q(rho*) < 0, everywhere on D**: proven exactly
  (Theorem 4.3, Corollary 4.4), via certified root isolation of a
  degree-1807 discriminant combined with a topological continuity
  argument.
- **Central theorem (V(z) strictly convex on D, for every t>1)**: proven,
  by the chain above (Sec. 5).

## References

Beyond the Hohmann-transfer context already cited in Theorems 1-2
(Prussing, 1992), this result uses only standard computer algebra
techniques (resultants, modular polynomial reduction, Bernstein-basis
positivity certificates, and certified ball-arithmetic root isolation),
not separately cited external theorems.

---

## Appendix A. The defining quartic F(t,w,rho), and an independent re-check of F4>0

```python
import sympy as sp
t, w, rho = sp.symbols('t w rho')

F4 = 3*t**4*(t - 1)**2*(t + 1)*(3*t**3*w**2 - 8*t**3*w + 6*t**3
     + 3*t**2*w**2 - 8*t**2*w + 6*t**2 + 3*t*w**2 - 8*t*w + 3*w**2)
F3 = 4*t**3*(t - 1)*(3*t**6*w**3 - 13*t**6*w**2 + 16*t**6*w - 6*t**6
     + 3*t**5*w**3 - 13*t**5*w**2 + 16*t**5*w - 6*t**5
     + 3*t**4*w**3 - 10*t**4*w**2 + 8*t**4*w - 6*t**4
     + 6*t**3*w**3 - 10*t**3*w**2 + 16*t**3*w - 6*t**3
     + 3*t**2*w**3 - 13*t**2*w**2 + 16*t**2*w
     + 3*t*w**3 - 13*t*w**2 + 3*w**3)
F2 = 2*t**2*(2*t**8*w**4 - 18*t**8*w**3 + 47*t**8*w**2 - 40*t**8*w + 6*t**8
     + 2*t**6*w**4 - 8*t**6*w**3 - 3*t**6*w**2 + 8*t**6*w + 4*t**6
     + 4*t**5*w**4 - 6*t**5*w**3 + 8*t**5*w**2 + 8*t**5*w
     - 6*t**4*w**3 - 3*t**4*w**2 + 6*t**4
     + 4*t**3*w**4 - 8*t**3*w**3 - 40*t**3*w
     + 2*t**2*w**4 + 47*t**2*w**2 - 18*t*w**3 + 2*w**4)
F1 = -4*t**2*(t - 1)*(2*t**2*w - 3*t**2 + 2*t*w - 3*t + 2*w)*(
     t**5*w**3 - 3*t**5*w**2 + 4*t**5*w - 2*t**5
     + t**3*w**3 - 2*t**3*w**2 - 2*t**3
     + t**2*w**3 + 4*t**2*w - 3*t*w**2 + w**3)
F0 = t**2*(t - 1)**2*(t**2*w**2 + 2*t**2 + w**2)*(
     2*t**2*w - 3*t**2 + 2*t*w - 3*t + 2*w)**2

F_expr = F4*rho**4 + F3*rho**3 + F2*rho**2 + F1*rho + F0
```

**Independent re-check.** F4 is a product of `3*t^4*(t-1)^2*(t+1)`
(strictly positive for t>1) and the bracket
B(t,w) = 3(t^3+t^2+t+1)*w^2 - 8t(t^2+t+1)*w + 6t^2(t+1)
(an equivalent regrouping, by w-degree, of the second factor displayed
above for F4). Individual monomials of B have mixed sign, so positivity
cannot be read off term by term; instead, B is a quadratic in w with
positive leading coefficient 3(t^3+t^2+t+1) for t>0, and discriminant
(8t(t^2+t+1))^2 - 4*3(t^3+t^2+t+1)*6t^2(t+1) = -8*t^2*(t-1)^2*(t^2+4t+1)
(verified symbolically), which is <=0 for every real t and strictly
negative whenever t != 1. Hence for every t>1, B(t,.) has no real root
and is strictly positive for every real w, so F4(t,w)>0 for every t>1,
w>0. A dense numerical sweep (t in [1.001,50], w in [0.001,50], 2000
random points) confirms F4>0 at every sampled point, with no exception,
consistent with this exact argument. Reconstructing RX, RY, RZ, N_rho's
full derivation, and the degree-1807 discriminant chain of Sec. 4, from
scratch has not been done (see the reproducibility note at the top of
this document).
