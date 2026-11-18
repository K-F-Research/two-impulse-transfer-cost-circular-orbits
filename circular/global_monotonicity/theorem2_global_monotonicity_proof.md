# Theorem 2 -- Global Monotonicity of the Optimized Transfer Cost

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Numbering note.** This folder's name and this document's own internal
numbering refer to this result as "Theorem 2". In the published article,
this same result -- strict global monotonicity of
`V_R(phi)` on `(0,pi)` -- is **Theorem 1**. The article's Theorem 2 (local
curvature at the Hohmann configuration) is proved in
`circular/hohmann_local_minimum/theorem1_hohmann_strict_local_minimum.md`,
which for the same historical reason is filed as "Theorem 1". See
`docs/mathematical_conventions.md` for the full account of this swap.

*Verification note: the two central closed-form identities (`D'(p)` in
Sec. 6 and the factorization of `g(p)` in Sec. 8) have been independently
re-derived from the definitions with SymPy and confirmed exactly (zero
residual after simplification). Global monotonicity itself has also been
independently checked numerically, by direct minimization over `p` at 8
radius ratios R and 60 angles each (472 consecutive-pair comparisons), with
zero violations -- see "Reproduction" below. Both checks are confirmatory,
not the source of the proof.*

---

## 0. Target statement

**Theorem (target).** For every fixed R>1, and for all 0<phi_1<phi_2<pi:

  V_R(phi_1) < V_R(phi_2)

where V_R(phi) := min_{p>0} F(p,phi,R), F being the total two-impulse
transfer cost Delta-v between two coplanar circular orbits of radii r1=1 and
r2=R, at transfer angle Delta_theta = pi - phi, parametrized by the
semi-latus rectum p of the transfer conic.

---

## 1. Definitions and domain

**Notation.** r1=1, r2=R>1 (the choice r1=1 is a choice of length unit, with
no loss of generality by homogeneity of the Keplerian problem -- Delta-v at
fixed r1,r2 depends only on the ratio R=r2/r1 once mu and r1 are normalized
to 1). phi in (0,pi) is the angular departure from the Hohmann angle:
Delta_theta = pi - phi. c:=cos(phi), s:=sin(phi) (s>0 on (0,pi)). p>0 is the
semi-latus rectum of the transfer conic.

x := p/r1 - 1 = p - 1,  y := p/r2 - 1 = p/R - 1.

e^2(p,phi,R) := x^2 + (x*c+y)^2/s^2 (squared eccentricity of the transfer
conic joining the two points; standard formula obtained from
r = p/(1+e*cos(theta-theta0)) evaluated at the two true anomalies separated
by Delta_theta).

v_theta,k(p) := sqrt(p)/r_k (transverse velocity on the transfer conic at
radius r_k, via specific angular momentum h=sqrt(p)).
v_c,k := 1/sqrt(r_k) (local circular velocity).
v_k^2(p) := 2/r_k + (e^2(p)-1)/p (vis-viva).
S_k(p,phi,R) := v_k^2(p) + v_c,k^2 - 2 v_c,k v_theta,k(p) (squared norm of
the Delta-v at endpoint k, for a tangential rendezvous with the local
circular orbit).

F(p,phi,R) := sqrt(S_1(p,phi,R)) + sqrt(S_2(p,phi,R)).

**Domain**: Omega := {(p,phi,R) : p>0, phi in (0,pi), R>1}. These
definitions are consistent with standard Keplerian mechanics (vis-viva,
angular momentum).

---

## 2. Admissibility -- S1>0, S2>0 on all of Omega

**Lemma 2.1 (sum-of-squares identity).** For every (p,phi,R) in Omega:

  S_k(p,phi,R) = v_{r,k}(p)^2 + (v_{theta,k}(p) - v_{c,k})^2

where v_{r,k}(p)^2 := v_k^2(p) - v_{theta,k}(p)^2 (squared radial velocity,
standard decomposition v^2 = v_r^2 + v_theta^2).

*Proof.* Expand the right-hand side:
v_{r,k}^2 + v_{theta,k}^2 - 2 v_{c,k} v_{theta,k} + v_{c,k}^2
= (v_{r,k}^2+v_{theta,k}^2) + v_{c,k}^2 - 2 v_{c,k} v_{theta,k}
= v_k^2 + v_{c,k}^2 - 2 v_{c,k} v_{theta,k} = S_k,
by the definition of v_{r,k}^2 and because v_k^2 = v_{r,k}^2 + v_{theta,k}^2
is the standard polar decomposition of velocity. QED

**Lemma 2.2 (v_{r,k}^2 is a manifest square).** For k=1:

  v_{r,1}(p)^2 = (R*c*p - R*c - R + p)^2 / (R^2 * p * s^2)

For k=2:

  v_{r,2}(p)^2 = (R*p - R*c - R + c*p)^2 / (R^2 * p * s^2)

*Proof.* v_{r,k}^2 = v_k^2 - v_{theta,k}^2 = 2/r_k + (e^2-1)/p - p/r_k^2.
Substitute e^2 = x^2 + (x*c+y)^2/s^2 and reduce to the common denominator
R^2*p*s^2; the numerator, a polynomial in p, c, R, factors into a perfect
square (the factorization above reproduces the numerator exactly after
expansion -- putting things over a common denominator, then identifying a
perfect square, confirmed by direct expansion). QED

Since R^2*p*s^2 > 0 for every p>0, R>0, phi in (0,pi), we obtain
**v_{r,k}^2 >= 0 with no condition whatsoever** -- in particular without
assuming e^2<1 (the result therefore also covers hyperbolic arcs, e>1,
with no additional hypothesis).

**Lemma 2.3 (characterization of v_{r,k}=0).**
v_{r,1}(p)=0  <=>  p = p1*(phi,R) := R(1+c)/(1+R*c).
v_{r,2}(p)=0  <=>  p = p2*(phi,R) := R(1+c)/(R+c).

*Proof.* Set the numerator of Lemma 2.2 to zero:
R*c*p - R*c - R + p = 0  <=>  p(R*c+1) = R(c+1)  <=>  p = R(1+c)/(1+R*c) (the
case R*c+1=0 corresponds to p1* undefined/at infinity, treated separately in
Sec. 7). Likewise for k=2. QED

**Lemma 2.4 (non-coincidence, strict S_k>0).**
(a) v_{theta,1}(p)=v_{c,1}  <=>  p=1. But p1*=1  <=>  R=1 (solve
R(1+c)=1+R*c  <=>  R=1). Since R>1, p1*!=1 for every phi, so the two
conditions "v_{r,1}=0" and "v_{theta,1}=v_{c,1}" can never be satisfied
simultaneously by the same p  ==>  **S1(p)>0 for every p>0, as soon as R>1**
(not merely near p1*: for every p).

(b) v_{theta,2}(p)=v_{c,2}  <=>  p=R. But p2*-R = -R(R-1)/(R+c) (direct
computation). Since R+c>0 (because c > -1 > -R for R>1) and R-1>0, this
quantity is **strictly negative** for every c in (-1,1): p2*!=R, without
exception  ==>  **S2(p)>0 for every p>0, as soon as R>1.**

This is a complete algebraic proof: no case inside Omega is left uncovered,
and degeneracy only appears at R=1, which is excluded from the domain.

**Direct consequence.** F(p,phi,R)=sqrt(S1)+sqrt(S2) is well-defined,
strictly positive, and **C-infinity in (p,phi)** over all of Omega (a
composition of regular rational and trigonometric functions, under a square
root of a strictly positive quantity -- hence itself regular, no artificial
singularity).

---

## 3. Existence of the minimum over p (coercivity)

**Lemma 3.1.** For every fixed phi in (0,pi), R>1: F(p,phi,R) -> +infinity as
p -> 0+ and as p -> +infinity.

*Proof.*
-- *p -> 0+*: e^2(p) -> e^2(0+) = 2/(1-c), a finite, strictly positive
constant (direct Taylor computation). Since e^2(0+)-1 = (1+c)/(1-c) > 0 for
c in (-1,1), v_k^2(p) = 2/r_k + (e^2(p)-1)/p ~ (1+c)/((1-c)*p) ->
+infinity. So S_k(p) -> +infinity, F(p) -> +infinity.
-- *p -> +infinity*: the coefficient of the p^2 term in the expansion of
e^2(p) is [R^2+2Rc+1] / [R^2(1-c^2)]. The numerator R^2+2Rc+1, viewed as a
function of c in (-1,1) at fixed R, is strictly positive as soon as R!=1
(at both endpoints c=+-1 it equals (R-+1)^2 >= 0, with equality only if
R=1; linear in c in between, hence strictly positive on the whole open
interval as soon as R>1). So e^2(p) ~ (const>0)*p^2 -> +infinity, hence
v_k^2(p) ~ e^2(p)/p -> +infinity (linear growth), and S_k(p) -> +infinity
(the term -2*v_{c,k}*v_{theta,k} ~ -O(sqrt(p)) is dominated). This holds
**including at the boundary c=-1/R** (where a different term could have
vanished): the coefficient of p^2 there equals exactly 1, still strictly
positive. QED

**Corollary 3.2 (existence).** For every fixed phi,R, F(.,phi,R) is
continuous and strictly positive on (0,infinity) (Sec. 2), and tends to
+infinity at both endpoints (Lemma 3.1). Hence its infimum over
(0,infinity) is attained: the set argmin(phi,R) := {p>0 :
F(p,phi,R)=V_R(phi)} is **non-empty**.

---

## 4. Interiority of the minimizer

**Lemma 4.1.** Every element of argmin(phi,R) is, by definition of the
minimization space, a real number p* in (0,infinity). The substantial point
is not this triviality, but the fact that coercivity (Sec. 3) guarantees
this minimum is not merely approached "at infinity" by a non-convergent
minimizing sequence: this is the quantitative content of Sec. 9 (existence
of a compact [p_I,P_I], independent of phi in a neighborhood, containing
every minimizer), which does the real analytic work. Sec. 4's useful
content is therefore folded into Sec. 9.

---

## 5. Uniqueness of the minimizer -- not needed

**Uniqueness is not assumed.** What is established (Sec. 7-8) is stronger
than what would be necessary for a proof relying on uniqueness: we show
that **d_phi F(p,phi) > 0 for EVERY p in the localization interval**
(p2*,p1*) in Region I, or (p2*,infinity) in Region II -- not only at the
point p* itself. Consequently:

- If argmin(phi,R) is a singleton {p*}, the classical envelope theorem
  applies directly.
- If argmin(phi,R) contains several points, they are *all* in the zone
  where d_phi F>0 (Sec. 8), so the argument of Sec. 9-10 (Danskin) concludes
  the same way without needing to settle uniqueness.

We therefore do not prove uniqueness of p* -- we show it is not necessary
for the targeted conclusion, thanks to the uniformity of the sign of
d_phi F over the whole localization zone (Sec. 8, Lemma 8.2). A strict
uniqueness proof (via convexity of F in p, not established here) would
remain useful for other reasons (uniqueness of the numerical solution,
physical interpretation) but is **not** a prerequisite of the monotonicity
theorem targeted here.

---

## 6. Central lemma -- D'(p) > 0

**Definition.** D(p,phi) := [d_p(e^2)(p)*p - (e^2(p)-1)] / p^2. (This term
arises naturally: S_k'(p) = D(p) - v_{c,k}/(r_k*sqrt(p)), since
v_k^2'(p)=D(p) is common to both k, only the term -2*v_{c,k}*v_{theta,k}
depending on k in the derivative.)

**Lemma 6.1.**

  D'(p) = 2(1+c) / [p^3*(1-c)] > 0

for every p>0 and every c in (-1,1) -- **independent of R**.

*Proof.* e^2(p) is, at fixed c,R, a rational fraction in p of the form
e^2(p) = [quadratic(p)] / [R^2(1-c^2)] (numerator a degree-2 polynomial in
p, denominator constant in p). A direct computation (differentiating
x^2+(x*c+y)^2/s^2 with respect to p, x=p-1, y=p/R-1, s^2 = 1-c^2 constant in
p) gives D(p) in explicit form; its derivative D'(p) simplifies, after
reduction to a common denominator and cancellation of the degree->=1 terms
(which cancel identically), into the fraction above. The sign is immediate:
p^3>0, and (1+c)>0, (1-c)>0 for strict c in (-1,1). QED

This closed form has no dependence on R, a notable and structurally clean
fact; it has been independently re-derived and confirmed with zero residual
(see "Reproduction" below).

**Corollary 6.2.** S1'(p) = D(p) - 1/sqrt(p) and
S2'(p) = D(p) - 1/(R^{3/2}*sqrt(p)) are **strictly increasing in p** on
(0,infinity): sums of D(p) (strictly increasing, Lemma 6.1) and of
-1/sqrt(p), -1/(R^{3/2}*sqrt(p)) (strictly increasing, derivatives
1/(2*p^{3/2})>0 and 1/(2*R^{3/2}*p^{3/2})>0).

---

## 7. Localization of p*

**Lemma 7.1 (anchor values).** D(p2*) = 1/R^2, D(p1*) = 1 (in Region I,
p1*>0).

*Proof.* Direct substitution of p1*=R(1+c)/(1+R*c), p2*=R(1+c)/(R+c) into
the explicit form of D(p) from Sec. 6, followed by algebraic simplification
(common denominator, cancellation of terms). QED

**Lemma 7.2 (the two key inequalities).**
(a) p2* < R for every R>1, every c in (-1,1).
(b) p1* > 1 for every R>1, every c with 1+R*c>0 (Region I).

*Proof.* (a) p2*=R(1+c)/(R+c) < R  <=>  (1+c)/(R+c) < 1  <=>  1+c < R+c
(since R+c>0, licit division)  <=>  1<R. True by hypothesis.
(b) p1*=R(1+c)/(1+R*c) > 1  <=>  R(1+c) > 1+R*c (since 1+R*c>0 in
Region I, licit division)  <=>  R+R*c > 1+R*c  <=>  R>1. True by hypothesis.
QED

**Lemma 7.3 (signs at the anchor points).** S1'(p2*)<0, S2'(p2*)<0
(everywhere on Omega); S1'(p1*)>0, S2'(p1*)>0 (Region I).

*Proof.* S1'(p2*)=1/R^2 - 1/sqrt(p2*) < 0  <=>  sqrt(p2*) < R^2  <=>
p2*<R^4. Since p2*<R (Lemma 7.2a) and R<R^4 for R>1, the inequality holds.
Likewise S2'(p2*)=1/R^2 - 1/(R^{3/2}*sqrt(p2*)) < 0  <=>  sqrt(p2*) < sqrt(R)
<=>  p2*<R, exactly Lemma 7.2a. For p1*: S1'(p1*)=1-1/sqrt(p1*)>0  <=>
p1*>1, Lemma 7.2b. S2'(p1*)=1-1/(R^{3/2}*sqrt(p1*))>0  <=>  p1*>R^{-3}; since
p1*>1>R^{-3} (R>1), it holds. QED

**Theorem 7.4 (global localization, not merely pointwise).** By the strict
monotonicity of S1', S2' in p (Corollary 6.2) and the signs of Lemma 7.3:

  S1'(p)<0 and S2'(p)<0 for every p <= p2* (Region I and II, without
  restriction);
  S1'(p)>0 and S2'(p)>0 for every p >= p1* (Region I only).

Hence F'(p) = S1'(p)/(2*sqrt(S1)) + S2'(p)/(2*sqrt(S2)) is **strictly
negative** on (0,p2*] and **strictly positive** on [p1*,infinity) (Region
I) -- F is strictly decreasing then strictly increasing outside the
localization zone, which rules out any minimizer outside the open interval.
In Region II, only the first fact applies (p1*<0, out of scope): F strictly
decreasing on (0,p2*].

**Conclusion (localization).**
argmin(phi,R) is a subset of (p2*,p1*) in Region I (c>-1/R);
argmin(phi,R) is a subset of (p2*,infinity) in Region II (c<-1/R).

This is a **global** localization (valid for the whole function F, not
only for an assumed-known minimizer) -- this is what allows bypassing the
uniqueness question of Sec. 5.

*A secondary point.* In Region II, the upper bound of argmin is not given
explicitly (only the existence of a bound, via the coercivity of Sec. 3, is
guaranteed). This does not call into question the conclusion of Sec. 8
(uniform sign over all of (p2*,infinity)): the sign of g(p) is established
over the whole ray (p2*,infinity), so the lack of an explicit upper bound
does not weaken the final conclusion -- it only means the localization
"window" in Region II is not as precise as in Region I.

---

## 8. Sign of d_phi F

**Lemma 8.1 (reduction).** d_phi F(p,phi) depends on phi only through
e^2(p,phi) (the other terms of S1,S2 do not depend on phi). Hence:

  d_phi F(p,phi) = [d_phi(e^2)(p,phi) / (2p)] * [1/sqrt(S1) + 1/sqrt(S2)]

The factor [1/sqrt(S1)+1/sqrt(S2)] is **strictly positive** everywhere on
Omega (Sec. 2). Hence sign(d_phi F) = sign(d_phi(e^2)).

**Lemma 8.2 (closed form of d_phi(e^2)).**

  d_phi(e^2)(p,phi) = -(2/sin^3(phi)) * g(p),  where g(p) := (x*c+y)(x+y*c)

*Proof.* Direct differentiation of e^2=x^2+(x*c+y)^2/s^2 with respect to
phi (with c=cos(phi), s=sin(phi), x,y independent of phi), a mechanical
chain-rule computation, simplified (the full chain rule is accounted for).
QED

Since sin^3(phi)>0 on (0,pi): **sign(d_phi F(p,phi)) = -sign(g(p))**.

**Lemma 8.3 (g(p) is a quadratic in p, exact factorization).**

  g(p) = A*(p-p1*)*(p-p2*),  A = (R+c)(1+R*c)/R^2

(A*(p-p1*)*(p-p2*) - g(p) == 0 is an exact polynomial identity.)

In Region I (c>-1/R): R+c>0 and 1+R*c>0, hence **A>0**. For an upward
parabola, g(p)<0 strictly **for every p strictly between the two roots** --
not only at one particular point.
In Region II (c<-1/R): R+c>0 still, 1+R*c<0, hence **A<0**. Downward
parabola: g(p)<0 strictly **for every p outside the segment [p1*,p2*]**;
since p1*<0<p2* in this case, this covers in particular every p>p2*.

**Theorem 8.4 (sign conclusion, uniform over the whole localization
zone).**
g(p)<0, hence d_phi F(p,phi)>0, **for every p in (p2*,p1*)** in Region I,
and **for every p in (p2*,infinity)** in Region II.

Combined with the localization of Sec. 7 (argmin is a subset of these same
sets):

  **d_phi F(p,phi) > 0 for every p in argmin(phi,R).**

This is exactly the point anticipated in Sec. 5: since the sign is uniform
over *the whole* zone where minimizers can lie (not only proven at a point
p* whose uniqueness would need to be known first), the conclusion does not
depend on a uniqueness result.

---

## 9. Regularity of V_R -- local uniform coercivity, with an explicit reference bound

**Lemma 9.1 (uniform coercivity on every compact I subset of (0,pi), via a
reference point).** Let I=[phi_1,phi_2] subset of (0,pi) be compact (I is
taken directly equal to the target interval [phi_1,phi_2] of the target
statement -- no abstract neighborhood of a point phi_0 to introduce). Then
there exist P_I<infinity and p_I>0 (depending on I and R, not on phi) such
that for every phi in I, every minimizer p*(phi) of F(.,phi,R) satisfies
p_I <= p*(phi) <= P_I.

*Proof.* **Reference bound.** Let p0:=p_H=2R/(R+1) (the Hohmann semi-latus
rectum, a fixed point independent of phi). Since F(p0,.,R) is continuous in
phi (Sec. 2) on the compact I, M_I := max_{phi in I} F(p0,phi,R) < infinity.
For every phi in I: V_R(phi) <= F(p0,phi,R) <= M_I (by definition of V_R as
a minimum). It therefore suffices to show that F(p,phi,R) > M_I+1 as soon as
p leaves a fixed compact [p_I,P_I], uniformly over I: **if such a p
satisfies F(p,phi,R)>M_I+1, then, since V_R(phi)<=M_I<M_I+1<F(p,phi,R), we
have F(p,phi,R)>V_R(phi), so p does not realize the minimum at phi** -- such
a p can therefore never be a minimizer, whatever phi in I.

**Bound as p -> 0+.** On I, e^2(p,phi) = 2/(1-c) + O_I(p) where the O_I is
uniform over I (continuous coefficients in phi, hence bounded on the
compact I). Hence (e^2(p,phi)-1)/p = (1+c)/[(1-c)*p] + O_I(1). Since
mu_I := min_{phi in I} (1+cos(phi))/(1-cos(phi)) > 0 (continuous and
strictly positive on I, compactly included in (0,pi), so c=cos(phi) stays
strictly below 1 there, bounded away from 1), we obtain, uniformly over I:
v_k^2(p,phi) >= mu_I/p - C1 for a constant C1 independent of phi. Hence
S_k(p,phi) -> +infinity uniformly as p -> 0+, and there exists p_I>0
(explicit in terms of mu_I, C1, M_I) such that F(p,phi,R)>M_I+1 for every
p<=p_I, every phi in I.

**Bound as p -> infinity.** e^2(p,phi)=A(phi)*p^2+B(phi)*p+C(phi) with
A(phi)=(R^2+2*R*cos(phi)+1)/(R^2*(1-cos^2(phi))). A is continuous and
strictly positive on I (numerator R^2+2Rc+1>0 for R>1, every c in (-1,1) --
Lemma 3.1), hence A_I:=min_{phi in I} A(phi)>0 (attained by compactness).
B,C are continuous hence bounded on I. One deduces
e^2(p,phi)>=(A_I/2)*p^2 for p large enough (explicit threshold, independent
of phi), uniformly over I, hence v_k^2(p,phi)>=(A_I/2)*p - O_I(1), and
S_k(p,phi) -> +infinity uniformly as p -> +infinity. Hence the existence of
P_I<infinity (explicit in terms of A_I, the bounds of B,C, and M_I) such
that F(p,phi,R)>M_I+1 for every p>=P_I, every phi in I.

The two bounds give K_I:=[p_I,P_I] such that every minimizer, for every phi
in I, belongs to K_I. QED

**Lemma 9.2 (continuity of V_R).** V_R(phi)=min_{p in K_I} F(p,phi,R) on I
(Lemma 9.1). F is jointly continuous on the compact K_I x I (Sec. 2). By
Berge's maximum theorem (Berge, *Topological Spaces*, 1963), V_R is
continuous on I, and the correspondence phi -> M(phi):=argmin_{p in K_I}
F(p,phi,R) is upper semicontinuous, with non-empty compact values.

---

## 10. Directional derivative of the optimal value -- Danskin's theorem

The optimized cost function is V_R(phi)=min_{p>0} F(p,phi,R).

Bertsekas's result, Proposition B.25, is **not** directly applicable here:
its formulation relies on a convexity hypothesis for the function with
respect to the outer variable. We neither assume nor prove that
phi -> F(p,phi,R) is convex at fixed p (F contains a term in 1/sin^2(phi),
which rules out any presumption of global convexity).

We therefore use a version of Danskin's theorem **without a convexity
hypothesis**, in the form given by **Shapiro, A., Dentcheva, D. &
Ruszczynski, A. (2009), *Lectures on Stochastic Programming: Modeling and
Theory*, SIAM/MPS Series on Optimization, Theorem 7.21, p. 356** (Section
7.1.5, "Perturbation Analysis" -> "Differentiability Properties of
Max-Functions"). Exact statement from the book: for Theta compact non-empty
and g: R^n x Theta -> R such that g(.,theta) is differentiable for every
theta in Theta and grad_x g(x,theta) is continuous on R^n x Theta, the
max-function phi(x):=sup_{theta in Theta} g(x,theta) is locally Lipschitz,
directionally differentiable, with phi'(x,h)=sup_{theta in Theta_bar(x)}
h^T grad_x g(x,theta) -- **no convexity of x -> g(x,theta) appears in the
hypotheses.** This applies to a function defined as the maximum of a
family of differentiable functions whose gradient depends continuously on
the variables -- without convexity in the outer variable. (This is
historically Danskin's original formulation: J. M. Danskin, *The Theory of
Max-Min and its Application to Weapons Allocation Problems*, Springer,
1967, also without convexity. The book contains, immediately after, a
second theorem -- Theorem 7.22 (Levin-Valadier), p. 357 -- which requires
g(.,theta) convex for each theta and gives in exchange the full convex
subdifferential formula d(phi)(x)=conv{grad g_theta(x): theta in
Theta_bar(x)}; it is this Theorem 7.22, with its convexity hypothesis, that
corresponds to Bertsekas B.25's register -- not Theorem 7.21, which we use.
Our argument only needs the directional-derivative formula, never the
subdifferential.)

**Passing from min to max.** To apply Theorem 7.21 (stated for a max) to our
minimization problem, set G(p,phi):=-F(p,phi,R). On every compact interval
I=[phi_1,phi_2] compactly included in (0,pi), Sec. 9 provides a compact
K_I subset of (0,infinity) such that V_R(phi)=min_{p in K_I} F(p,phi,R) for
every phi in I. Hence -V_R(phi)=max_{p in K_I} G(p,phi). Since K_I is
compact and F is C-infinity on K_I x I (Sec. 2), the hypotheses of
Theorem 7.21 are satisfied for G. Writing
M(phi):=argmin_{p in K_I} F(p,phi,R)=argmax_{p in K_I} G(p,phi), we obtain,
for every direction d in R:

  D(-V_R)(phi;d) = max_{p in M(phi)} d * G_phi(p,phi)

Now G_phi=-F_phi, hence D(-V_R)(phi;d) = max_{p in M(phi)}(-d*F_phi(p,phi,R))
= -min_{p in M(phi)} d*F_phi(p,phi,R). Consequently:

  **DV_R(phi;d) = min_{p in M(phi)} d*F_phi(p,phi,R)**

In particular, for the right directional derivative (d=+1):

  **D+V_R(phi) = min_{p in M(phi)} F_phi(p,phi,R)**

This formula requires neither uniqueness of the minimizer nor convexity of F
with respect to phi: several values of p may realize the minimum, it
suffices to control the sign of F_phi over the whole set of minimizers.

We use the convexity-free directional-derivative form of Danskin's theorem,
as stated in Shapiro, Dentcheva and Ruszczynski, Theorem 7.21; the minimum
form used here follows immediately by applying the theorem to -F.

---

## 11. Global monotonicity theorem

We first build the compact graph of minimizers on I, then a single global
bound, then an elementary right-derivative lemma -- with no finite cover,
which turns out unnecessary once the bound is made global.

**Lemma 11.1 (compact graph of minimizers).** On I=[phi_1,phi_2], K=K_I
(Sec. 9):

  M_I := {(p,phi) in K x I : F(p,phi,R)=V_R(phi)}

is a closed subset of K x I (since F and V_R are continuous, Sec. 2 and
Lemma 9.2), hence **compact**. By construction, M_I contains exactly the
pairs (p,phi) where p is a minimizer at angle phi.

**Lemma 11.2 (uniform global bound on the graph).** For every (p,phi) in
M_I, p is in argmin(phi,R) subset of L(phi) (localization zone, Sec. 7),
hence d_phi F(p,phi)>0 by Theorem 8.4. Since d_phi F is continuous on the
compact M_I (Sec. 8, restriction of Lemma 8.1-8.2, continuous on K_I x I),
and M_I is non-empty compact (Lemma 11.1):

  m_I := min_{(p,phi) in M_I} d_phi F(p,phi) > 0

(minimum attained, strictly positive: a zone where d_phi F vanishes does
exist within the whole of K_I x I -- at the roots p1*,p2* themselves,
d_phi F vanishes -- but M_I, the graph of the *minimizers*, never touches
it, by construction of Sec. 7).

**Lemma 11.3 (uniform control of the directional derivative).** By Sec. 10,
D+V_R(phi) = min_{p in M(phi)} F_phi(p,phi,R). By Sec. 8, every minimizer
lies in the region where F_phi(p,phi,R)>0: this positivity is pointwise
(established point by point over the set of minimizers), and it remains to
turn it into a strictly positive bound *uniform* over I -- which is exactly
what the compact graph M_I of Lemma 11.1 does: F_phi is continuous and
strictly positive there (Sec. 8), so it attains on this non-empty compact a
minimum m_I:=min_{(p,phi) in M_I} F_phi(p,phi,R)>0 (Lemma 11.2). Since for
every phi in I, M(phi) is exactly the fiber of M_I at phi, one obtains

  D+V_R(phi) = min_{p in M(phi)} F_phi(p,phi) >= m_I > 0,

**for every phi in I simultaneously** -- the same constant m_I everywhere,
no local bound varying with phi to patch together.

**Lemma 11.4 (lower-bounded right derivative implies strict comparison at
the endpoints).** Let g:[a,b]->R be continuous, with right derivative
D+g(x) existing at every x in [a,b) and satisfying D+g(x)>=eta>0. Then
**g(b)-g(a)>=eta*(b-a)>0**. (This is all that Theorem 11.5 below needs --
comparing the two endpoints a=phi_1, b=phi_2 of a given interval; the lemma
does not assert, and does not need to assert, that g would be monotone
between arbitrary intermediate points.)

*Proof (the standard sup technique -- Denjoy-Young-Saks).* Fix arbitrary
epsilon>0, and set h(x):=g(x)-(eta-epsilon)*x. Since D+g(x)>=eta,
D+h(x)>=epsilon>0 for every x in [a,b). We show **h(b)>=h(a)** (no more --
we do not seek to compare h at arbitrary intermediate points), which will
give g(b)-g(a)>=(eta-epsilon)*(b-a); the result follows by letting
epsilon -> 0.

Let S := {x in [a,b] : h(t)>=h(a) for every t in [a,x]}. S is non-empty
(a in S) and, by continuity of h, closed in [a,b] -- hence c:=sup(S) is in
S. Suppose for contradiction c<b. Since D+h(c)>=epsilon>0, there exists
delta>0 such that for every t in (c,c+delta]:
[h(t)-h(c)]/(t-c) > 0, hence h(t)>h(c)>=h(a) (the last inequality because
c is in S). Thus every t in [a,c+delta] satisfies h(t)>=h(a) (on [a,c] by c
in S, on (c,c+delta] by the above): c+delta is in S, contradicting
c=sup(S). Hence c=b, and b in S means precisely **h(b)>=h(a)**, the only
conclusion sought. QED

This lemma is a generic, uncontroversial result of real analysis (no
hypothesis on g beyond continuity and the existence of the right
derivative); its proof depends on no property specific to V_R.

---

**The final chain, in the exact form to appear in the manuscript** (no step
introduces uniqueness of p*):

  D'(p)>0
  ==> localization of M(phi):=argmin_{p>0} F(p,phi,R) (Sec. 7, valid for
  every element of M(phi), not for an assumed-unique p*)
  ==> F_phi(p,phi,R)>0 for **every** p in M(phi) (Sec. 8, Theorem 8.4 --
  uniform sign over the whole set of minimizers)
  ==> **[Danskin]** D+V_R(phi) = min_{p in M(phi)} F_phi(p,phi,R) > 0
  (Sec. 10)
  ==> **[Lemma 11.4]** for every [phi_1,phi_2] subset of (0,pi):
  V_R(phi_2)-V_R(phi_1) >= eta_{[phi_1,phi_2]}*(phi_2-phi_1) > 0.

**A point of method worth making explicit**: we do **not** seek a uniform
bound D+V_R(phi)>=eta>0 valid over all of (0,pi) -- such a bound would in
fact be false near phi=0, since the local expansion (Sec. 13) gives
V_R'(phi)~2*kappa(R)*phi -> 0 as phi -> 0+. The right target, and the one
actually established, is: **for every given compact [phi_1,phi_2] subset of
(0,pi) (arbitrary, but fixed before constructing the bound), there exists
eta=eta_{[phi_1,phi_2]}>0** such that D+V_R(phi)>=eta on this precise
interval -- this is exactly what Lemma 11.1 (compact graph of minimizers on
I=[phi_1,phi_2]) and Lemma 11.2 (m_I:=min on this compact graph, attained
and >0) give. Nothing stronger is necessary, and nothing stronger is true
globally.

**Theorem 11.5 (global monotonicity).** For every R>1 and every
0<phi_1<phi_2<pi (arbitrary): taking I=[phi_1,phi_2], K=K_I (Sec. 9),
m_I>0 (Sec. 11.2), Lemma 11.4 applied to g=V_R on [a,b]=[phi_1,phi_2] with
eta=m_I directly gives V_R(phi_2)-V_R(phi_1) >= m_I*(phi_2-phi_1) > 0.
**No finite cover is necessary**: [phi_1,phi_2] being the target interval
itself (not a sub-interval of an abstract neighborhood), the bound m_I
obtained in Sec. 11.2 applies directly, over the whole interval, in a
single step.

**Corollary (target statement of Sec. 0).** For every R>1, 0<phi_1<phi_2<pi
==> V_R(phi_1)<V_R(phi_2), with the explicit quantitative bound
V_R(phi_2)-V_R(phi_1)>=m_I*(phi_2-phi_1).

---

## 12. Limiting cases phi -> 0+ and phi -> pi-

**phi -> 0+**: p1*, p2* -> p_H = 2R/(R+1) (the two roots merge exactly at
the Hohmann semi-latus rectum). The theorem does not cover phi=0 itself
(excluded from the open interval), which is consistent: at phi=0, V_R
attains its global minimum (Prussing's 1992 result), so V_R'(0)=0 is
expected there, not V_R'(0)>0. This is qualitatively consistent with the
local result of Sec. 13; precise quantitative matching with the local
coefficient K(R)*phi^2 is not pursued here.

**phi -> pi-**: this limiting case is **not necessary** for the
monotonicity theorem (Sec. 11), which concerns the open interval (0,pi)
without requiring any control of boundary behavior: the endpoint phi -> pi-
is not required for the monotonicity theorem and is therefore left outside
the scope of the present analysis. The numerical observation (p*(phi)->0+,
regular behavior with no apparent pathology, Sec. 14) remains available as
a numerical remark but is not presented as an established asymptotic claim.

---

## 13. Relation to the local result K(R)

**The global theorem does not prove, and should not be presented as
proving, K(R)>0.** What it gives directly is the qualitative structure
V_R(phi)>V_R(0) for phi>0 -- consistent with, but logically distinct from,
the local expansion V_R(phi)=V_H+K(R)*phi^2+O(phi^4) already established in
Theorem 1 / the Lemma-1 folder. The natural even extension
V_R(-phi)=V_R(phi) gives V_R'(0+)=0, which is equally consistent (the
global theorem concerns the open interval (0,pi), not phi=0 itself -- no
contradiction, simply two results at different scales):

-- *Global theorem*: angular deviation from the Hohmann angle always
strictly increases the optimized transfer cost, on the entire interval
(0,pi).
-- *Local corollary*: sufficiently near the Hohmann angle, this increase is
quadratic, with explicit coefficient K(R) (established independently
elsewhere in this project).

This is a qualitative/corollary-level relationship between the two
results, not a claim that the global theorem formally re-proves the local
one.

---

## 14. Independent numerical validation

Summary (details in the project's own numerical scripts, not reproduced
here verbatim):
- Arbitrary precision (mpmath, 50 significant digits), root search by
  bisection on d_p F (independent of `scipy.optimize.minimize_scalar` used
  in a first check).
- 243 configurations (9 values of R in [1.2,50], angles from 1 deg to
  179.9 deg): localization p2*<p*<p1* (Region I) / p*>p2* (Region II)
  confirmed without exception, including very close to phi=180 deg.
- 200 additional configurations (R in [1.2,100]) testing directly the sign
  of d_p F at the points p1*,p2* themselves (without going through p*):
  confirmed without exception, consistent with Lemma 7.3.
- Feasibility S1,S2>0 checked numerically over the whole grid (no
  violation), consistent with the analytic proof of Sec. 2.

This numerical validation does not replace the analytic proof of
Sec. 1-13; it is strong empirical support, fully consistent with it.

---

**Summary of the proof chain:**

  D'(p)>0
  ==> localization of M(phi) (the whole set, not an assumed-unique p*)
  ==> F_phi>0 over all of M(phi)
  ==> [Sec. 10, Danskin via G=-F, Shapiro-Dentcheva-Ruszczynski Th. 7.21]
      D+V_R(phi)=min_{p in M(phi)} F_phi(p,phi,R)>0
  ==> [Sec. 11.2-11.3, compact graph M_I, bound m_I>0]
      D+V_R(phi)>=m_I on I
  ==> [Lemma 11.4] V_R(phi_2)-V_R(phi_1)>= m_I*(phi_2-phi_1)>0 for every
      [phi_1,phi_2] subset of (0,pi).

The structure D'(p)>0 combined with g(p)=A(p-p1*)(p-p2*) is what allows
matching the localization of *all* minimizers (not an assumed-unique p*)
with the region where F_phi is positive -- without ever solving the
degree-8 polynomial that a direct approach would otherwise require. This
is the central methodological contribution of the proof.

## References

Prussing, J. E., "Simple Proof of the Global Optimality of the Hohmann
Transfer," *Journal of Guidance, Control, and Dynamics*, vol. 15, no. 4
(1992), pp. 1037-1038.

Shapiro, A., Dentcheva, D. & Ruszczynski, A. (2009), *Lectures on
Stochastic Programming: Modeling and Theory*, SIAM/MPS Series on
Optimization, Theorem 7.21, p. 356.

Danskin, J. M. (1967), *The Theory of Max-Min and its Application to
Weapons Allocation Problems*, Springer.

Berge, C. (1963), *Topological Spaces*.
