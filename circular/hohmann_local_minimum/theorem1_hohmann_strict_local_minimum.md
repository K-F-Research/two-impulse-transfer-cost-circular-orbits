# The Hohmann Transfer as a Strict Local Minimum: Formal Proof

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Numbering note.** This folder's name and this document's own internal
numbering refer to this result as "Theorem 1". In the published article,
this same result -- strict positivity of the local
curvature coefficient at the Hohmann configuration -- is **Theorem 2**.
The article's Theorem 1 (global monotonicity of `V_R(phi)`) is proved in
`circular/global_monotonicity/theorem2_global_monotonicity_proof.md`,
which for the same historical reason is filed as "Theorem 2". See
`docs/mathematical_conventions.md` for the full account of this swap.

This is a complete, hand-verifiable proof: definitions, numbered statements,
and full proofs checkable by hand throughout. Confirmatory symbolic/numerical
checks that were run alongside this derivation are noted as such at the end
of each section; they are confirmation only, never the source of the results
stated here.

The notational conventions used throughout (`phi := pi - Delta_theta`, and
`R := r2/r1` in Section 8) are consistent with
`docs/mathematical_conventions.md`.

---

## 0. Scope of this document

This document establishes, by a complete and self-contained mathematical
argument (no step relies on computer-assisted numerical or symbolic
computation, except for a few final confirmations explicitly flagged as
such), the following result:

> **Theorem.** The Hohmann transfer between two coplanar circular orbits of
> radii $r_1 \neq r_2$ is not only the global minimum of the cost $\Delta v_T$
> among all two-impulse transfers (Prussing, 1992): it is also a **strict**
> local minimum in the sense of the transfer angle $\Delta\theta$ -- the cost
> increases strictly as soon as one departs, even slightly, from
> $\Delta\theta = \pi$.

The proof decomposes into five numbered lemmas (Lemma 1 to Lemma 5) plus one
local characterization ("Lemma 3" in what follows, for ease of reference --
not an independent additional statement, but the name given to the
characterization of $p^*(\Delta\theta)$ near $\pi$):

- **Lemma 1** (Sec. 2) and **Lemma 2** (Sec. 3).
- **"Lemma 3"** (characterization of $p^*$ near $\pi$, Sec. 4-8), including
  the justification, in Sec. 4bis, that the family $p-p_H=O(\varphi^2)$ is
  the only competitive family (not an unproven ansatz).
- **Lemmas 4-5** (symmetry, scale invariance, Sec. 8).
- **Theorem** (Sec. 9): conclusion.

## 1. Notation and framework

We adopt Prussing's (1992) formalism, to remain directly comparable with the
reference on which this work builds.

**Definition 1 (transfer orbit).** Let $r_1, r_2 > 0$, $r_1 \neq r_2$, be the
radii of two coplanar circular orbits around a central body of gravitational
parameter $\mu > 0$. A *transfer orbit* between them is a conic with focus at
the center, of semi-latus rectum $p>0$ and eccentricity $e \geq 0$, passing
through a point at distance $r_1$ (departure) and a point at distance $r_2$
(arrival), separated by a transfer angle $\Delta\theta \in (0,\pi]$ measured
along the direction of travel. Let $\nu_k$ denote the true anomaly at the
$k$-th passage point ($k=1,2$), so that the orbit equation gives

$$r_k = \frac{p}{1+e\cos\nu_k} \quad\Longleftrightarrow\quad e\cos\nu_k = \frac{p}{r_k}-1 =: \xi_k .$$

Note that $\xi_k$ depends only on $(p,r_k)$ -- not on $e$ alone, nor on
$\Delta\theta$.

**Definition 2 (velocities, transfer cost).** On the transfer orbit, the
radial and transverse velocity components at the $k$-th point are

$$v_{r,k} = \sqrt{\frac{\mu}{p}}\, e\sin\nu_k, \qquad v_{\theta,k} = \sqrt{\frac{\mu}{p}}\,(1+e\cos\nu_k) = \frac{\sqrt{\mu p}}{r_k}.$$

The local circular velocity is $v_{c,k} = \sqrt{\mu/r_k}$, purely tangential
(circular orbit). The impulse required at the $k$-th point is simply the
vector difference between the circular velocity $\mathbf v_{c,k}=(0,\,v_{c,k})$
(radial, transverse coordinates) and the velocity on the transfer orbit
$\mathbf v_k=(v_{r,k},\,v_{\theta,k})$ at that same point -- no choice or
optimization of impulse direction is involved, only the kinematics of the
instantaneous matching of the two velocities at a given point. Its norm is
given by the law of cosines (Pythagoras applied to $\mathbf v_k-\mathbf
v_{c,k}$):

$$\Delta v_k^2 = v_k^2 + v_{c,k}^2 - 2\, v_{c,k}\, v_{\theta,k}, \qquad v_k^2 := v_{r,k}^2+v_{\theta,k}^2 .$$

The total transfer cost is $\Delta v_T(p,\Delta\theta) := \Delta v_1 + \Delta v_2$
(the parameter $p$ and the angle $\Delta\theta$ determine the transfer orbit
up to a reflection, unambiguously for $\Delta\theta \in (0,\pi)$).

**Remark (fact established elsewhere, cited, not re-proven here).** Expanding
$v_k^2$ using Definition 1,

$$v_k^2 = v_{r,k}^2 + v_{\theta,k}^2 = \frac{\mu}{p}\Big(e^2\sin^2\nu_k + (1+e\cos\nu_k)^2\Big) = \frac{\mu}{p}\big(1+e^2+2e\cos\nu_k\big) = \frac{\mu}{p}\big(e^2 + 1+2\xi_k\big),$$

substituting $e\cos\nu_k=\xi_k$. Since $\xi_k = p/r_k - 1$ depends only on
$(p,r_k)$, **we obtain**

$$v_k^2(p,e^2) = \frac{\mu}{p}\,e^2 \;+\; \underbrace{\frac{\mu}{p}(1+2\xi_k)}_{\text{independent of } e},$$

that is, $v_k^2$, at fixed $p$, is an **affine** function of $e^2$, with slope
$\mu/p$. This is the vis-viva identity rewritten in these variables; it is
immediate and needs no separate proof.

## 2. Lemma 1: strict monotonicity of the cost in $e^2$, at fixed $p$

> **Lemma 1.** For fixed $p>0$ and $r_1,r_2>0$, the function
> $e^2 \mapsto \Delta v_T(p,e^2)$ is strictly increasing:
> $$\frac{\partial \Delta v_T}{\partial (e^2)}\bigg|_p \;=\; \frac{\mu}{2p}\left(\frac{1}{\Delta v_1}+\frac{1}{\Delta v_2}\right) \;>\; 0,$$
> provided $\Delta v_1>0$ and $\Delta v_2>0$ (non-degenerate transfer, i.e. the
> transfer orbit does not coincide with either circular orbit).

**Proof.** By Definition 2 and the remark in Sec. 1,

$$\Delta v_k^2 = v_k^2 + v_{c,k}^2 - 2v_{c,k}v_{\theta,k}.$$

Now $v_{c,k}=\sqrt{\mu/r_k}$ and $v_{\theta,k}=\sqrt{\mu p}/r_k$ depend neither
on $e$ nor on $e^2$ -- only on $p$ and $r_k$, both held fixed here. So the
only term on the right-hand side that depends on $e^2$ is $v_k^2$, and by the
remark in Sec. 1 this dependence is affine, with slope $\mu/p$:

$$\frac{\partial (\Delta v_k^2)}{\partial (e^2)}\bigg|_p = \frac{\partial v_k^2}{\partial (e^2)}\bigg|_p = \frac{\mu}{p}.$$

Since $\Delta v_k > 0$ by hypothesis, we may differentiate $\Delta v_k =
\sqrt{\Delta v_k^2}$:

$$\frac{\partial \Delta v_k}{\partial (e^2)}\bigg|_p = \frac{1}{2\Delta v_k}\cdot\frac{\partial(\Delta v_k^2)}{\partial(e^2)}\bigg|_p = \frac{\mu}{2p\,\Delta v_k} \;>\;0,$$

the sign being strict since $\mu>0$, $p>0$, $\Delta v_k>0$. Summing over
$k=1,2$ (linearity of the derivative, $\Delta v_T=\Delta v_1+\Delta v_2$):

$$\frac{\partial \Delta v_T}{\partial(e^2)}\bigg|_p = \frac{\mu}{2p\,\Delta v_1} + \frac{\mu}{2p\,\Delta v_2} = \frac{\mu}{2p}\left(\frac{1}{\Delta v_1}+\frac{1}{\Delta v_2}\right) > 0. \qquad\blacksquare$$

**Scope and hypotheses.** The hypothesis $\Delta v_k>0$ excludes only the
degenerate case where a point of the transfer orbit coincides in velocity
with the local circular orbit (zero impulse at that point) -- a measure-zero
set in parameter space, with no bearing on what follows. No other hypothesis
is used: the result holds for every $p>0$, every $r_1\neq r_2>0$, with no
restriction on $\Delta\theta$.

*Confirmatory check (not constitutive of the proof).* The formula in Lemma 1
was tested by symbolic computation and cross-checked against an independent
numerical evaluation on several cases. This test confirms the statement
above; it does not substitute for it.

## 3. Lemma 2: reduction of the sign of $\partial\Delta v_T/\partial\Delta\theta$ to that of $\partial(e^2)/\partial\Delta\theta$

> **Lemma 2.** At fixed $p>0$, for every $\Delta\theta$ where
> $e^2(p,\Delta\theta)$ is differentiable and the transfer is non-degenerate
> in the sense of Lemma 1,
> $$\operatorname{sign}\left(\frac{\partial \Delta v_T}{\partial \Delta\theta}\bigg|_p\right) = \operatorname{sign}\left(\frac{\partial (e^2)}{\partial \Delta\theta}\bigg|_p\right).$$

**Proof.** At fixed $p$, $\Delta v_T$ depends on $\Delta\theta$ only through
$e^2(p,\Delta\theta)$: indeed $v_{c,k}$ and $v_{\theta,k}$ do not depend on
$\Delta\theta$ at fixed $p,r_k$ (Definition 2), so, by the proof of Lemma 1,
$\Delta v_k$ depends on $\Delta\theta$ (at fixed $p$) only through $e^2$. The
chain rule then gives

$$\frac{\partial \Delta v_T}{\partial \Delta\theta}\bigg|_p = \frac{\partial \Delta v_T}{\partial (e^2)}\bigg|_p \cdot \frac{\partial (e^2)}{\partial \Delta\theta}\bigg|_p.$$

By Lemma 1, the first factor is strictly positive. The sign of the product is
therefore exactly that of the second factor. $\blacksquare$

**Consequence.** The whole question of the monotonicity of the cost in
$\Delta\theta$ (at fixed $p$, then, by the envelope theorem, in a
neighborhood of every point where $p^*(\Delta\theta)$ is differentiable)
reduces to the study of the sign of a single quantity:
$\partial(e^2)/\partial\Delta\theta$. This is a reduction result, not yet the
final result -- this sign still needs to be addressed, which is the object
of Lemma 3.

## 4. Lemma 3 (precise statement)

Lemma 2 reduces everything to the study of $e^2(p,\Delta\theta)$. By
Definition 1 applied at the two passage points, with $\xi_k = p/r_k-1$ and
writing $x:=\xi_1$, $y:=\xi_2$, eliminating $\nu_1,\nu_2$ between the two
relations $e\cos\nu_1=x$, $e\cos\nu_2=y$, $\nu_2-\nu_1=\Delta\theta$ gives the
identity, valid for every $\Delta\theta\notin\{0,\pi\}$:

$$e^2\sin^2\Delta\theta \;=\; x^2+y^2-2xy\cos\Delta\theta. \tag{$\ast$}$$

*(This identity -- a "law of cosines" form -- is equivalent to the fractional
formula $e^2 = x^2+\left(\frac{x\cos\Delta\theta-y}{\sin\Delta\theta}\right)^2$;
expanding the square and using $\sin^2\Delta\theta=1-\cos^2\Delta\theta$
passes from one to the other -- a direct algebraic check, not reproduced
here for brevity.)*

The right-hand side of $(\ast)$ is regular (no denominator) at
$\Delta\theta=\pi$, unlike the fractional form -- which is precisely why this
identity is the chosen starting point for what follows: it separates the
question of sign from the purely technical one of the apparent singularity
at $\sin\Delta\theta=0$.

What remains to be proven, to close the local argument announced in Sec. 0,
is the following statement.

> **Lemma 3 (to be proven).** There exists $\delta>0$ such that, for
> $\Delta\theta \in (\pi-\delta,\pi)\cup(\pi,\pi+\delta)$ (outside the
> physical domain $(0,\pi]$ for the second half, to be handled by symmetry),
> at the point $p=p^*(\Delta\theta)$ realizing $\min_p \Delta v_T(p,\Delta\theta)$,
> $$\frac{\partial (e^2)}{\partial \Delta\theta}\bigg|_{p=p^*(\Delta\theta)} \cdot \operatorname{sign}(\pi-\Delta\theta) \;>\;0,$$
> that is, $e^2$, evaluated along the optimal path $p^*(\Delta\theta)$,
> decreases strictly as $\Delta\theta \to \pi$.

This is the statement which, combined with Lemmas 1 and 2, would directly
give the result targeted in Sec. 0. It is also the real difficulty of the
problem: $e^2(p,\Delta\theta)$ at **fixed** $p$ (not at the optimal point)
has no uniform finite limit as $\Delta\theta\to\pi$ except exactly at
$p=p_{\text{Hohmann}}=2r_1r_2/(r_1+r_2)$ -- the optimal trajectory
$p^*(\Delta\theta)$ therefore crosses a region where a naive term-by-term
differentiation argument does not apply without further justification.

**Chosen proof strategy.** Rather than attempting to differentiate
$e^2(p^*(\Delta\theta),\Delta\theta)$ directly at $\Delta\theta=\pi$, treat
the Hohmann point as a fold-type singularity: set $\varphi=\pi-\Delta\theta$
and seek a local expansion $p^*(\varphi) = p_{\text{Hohmann}} + \varphi^2
q^\ast + O(\varphi^4)$. Sec. 4bis below establishes that this form is not an
ansatz adopted for convenience, but the only family of paths $p(\varphi)$
able to compete for the minimum; Sec. 5-8 then derive, entirely by hand, the
explicit value of $q^\ast$ and the conclusion.

## 4bis. The family $p-p_H=O(\varphi^2)$ is the only competitive one

Before restricting the search for the minimum to paths of the form
$p(\varphi)=p_H+\varphi^2q+O(\varphi^4)$, one must rule out that a path of a
different order -- in particular of order $\varphi^1$ -- could give a lower
cost. This closes in a few lines.

**Proposition (non-competitiveness of order-$\varphi$ paths).** Let
$p(\varphi)=p_H+\alpha\varphi+O(\varphi^2)$ with $\alpha\neq0$ be a family of
paths passing through $p_H$ at $\varphi=0$. Then
$\Delta v_T(p(\varphi),\pi-\varphi) \to \Delta v_{\text{Hohmann}}$ is **not**
attained at order $\varphi^2$: the gap to the Hohmann minimum along such a
path is of order $O(1)$ in $e^2$, hence strictly larger, for $\varphi$ small
enough, than any value reached by the family $p-p_H=O(\varphi^2)$.

*Proof.* With $\delta:=p-p_H=\alpha\varphi+O(\varphi^2)$, the formulas of
Sec. 5 for $x+y$ and $xy$ remain valid (they only assume $x,y$ affine in $p$,
hence in $\delta$, whatever the order of $\delta$ in $\varphi$):
$x+y=\delta s$, $xy=-e_H^2+O(\delta)$. First rewrite the right-hand side of
$(\ast)$ as $x^2+y^2-2xy\cos\Delta\theta = (x+y)^2-2xy(1+\cos\Delta\theta)$
(immediate algebraic identity). With
$1+\cos\Delta\theta=1-\cos\varphi=\varphi^2/2+O(\varphi^4)$ and
$xy=-e_H^2+O(\delta)$:
$$x^2+y^2-2xy\cos\Delta\theta = \delta^2s^2 \;-\;2\big(-e_H^2+O(\delta)\big)\Big(\frac{\varphi^2}2+O(\varphi^4)\Big) = \delta^2s^2 + e_H^2\varphi^2 + O(\delta\varphi^2,\varphi^4).$$
With $\delta=\alpha\varphi+O(\varphi^2)$ (so $\delta^2=\alpha^2\varphi^2+O(\varphi^3)$,
$\delta\varphi^2=O(\varphi^3)$) and $\sin^2\Delta\theta=\varphi^2+O(\varphi^4)$
on the left-hand side of $(\ast)$, identification at order $\varphi^2$ gives
$e^2\varphi^2=(\alpha^2s^2+e_H^2)\varphi^2+O(\varphi^3)$, i.e.
$$e^2(\varphi) \;\longrightarrow\; e_H^2+\alpha^2s^2 \;>\;e_H^2 \qquad(\varphi\to0),$$
a limit **strictly above** $e_H^2$, and nonzero as soon as $\alpha\neq0$ --
whereas, on the family $\delta=O(\varphi^2)$ (Sec. 5), $e^2\to e_H^2$ exactly,
with a gap of order $\varphi^2$ only (not a shifted limit). By Lemma 1
($\Delta v_T$ strictly increasing in $e^2$ at fixed $p$, hence a fortiori
increasing with the limiting value of $e^2$ along a path), the cost along an
order-$\varphi^1$ path tends, as $\varphi\to0$, to a value strictly above
$\Delta v_{\text{Hohmann}}$, with a gap that does not vanish with $\varphi$
-- hence, for $\varphi$ small enough, strictly worse than any path in the
family $\delta=O(\varphi^2)$, whose gap to $\Delta v_{\text{Hohmann}}$ tends
to $0$. $\blacksquare$

**Parity.** $\Delta v_T(p,\pi-\varphi)$ and $\Delta v_T(p,\pi+\varphi)$
correspond to transfers that are mirror images of one another (reversing the
direction of travel of the angle $\Delta\theta$ around $\pi$, at fixed $p$):
the function $\varphi\mapsto\min_p\Delta v_T(p,\pi-\varphi)$ is therefore
even. This rules out, on the competitive family, any term of order
$\varphi^1,\varphi^3,\ldots$ in the expansion of $p^*(\varphi)-p_H$
compatible with a regular minimum: the first available nontrivial order is
$\varphi^2$.

**Conclusion.** The search for the minimum can be restricted, without loss of
generality, to the family $p(\varphi)=p_H+\varphi^2q+O(\varphi^4)$: this is
no longer an ansatz, it is a consequence of Lemma 1 and parity. Sections
5-8 now determine $q^\ast$ and verify its positivity within this family,
whose potential optimality relative to any other family is already secured.

## 5. Expansion of $e^2$ near the Hohmann point

Set $\varphi=\pi-\Delta\theta$ and, following the announced ansatz,
$p=p_{\text{Hohmann}}+\delta$ with $\delta=\varphi^2 q$ for a free parameter
$q\in\mathbb R$ to be determined. Write $s:=\dfrac1{r_1}+\dfrac1{r_2}>0$; by
definition $p_{\text{Hohmann}}\,s = 2$.

**Step 1 -- $x+y$ and $xy$ in terms of $\delta$.** By Definition 1,
$x=p/r_1-1$ and $y=p/r_2-1$ are affine in $p$, hence in $\delta$:
$$x=x_H+\delta/r_1,\qquad y=y_H+\delta/r_2,\qquad x_H=e_H,\ y_H=-e_H,\quad e_H:=\frac{r_2-r_1}{r_1+r_2}.$$
($x_H,y_H$ are the values at $\delta=0$, computed directly from
$p_{\text{Hohmann}}/r_1-1$ and $p_{\text{Hohmann}}/r_2-1$.) From this, by a
direct computation (using $p_{\text{Hohmann}}s=2$):
$$x+y = \delta\, s, \qquad xy = -e_H^2 - \delta\, e_H^2 s + O(\delta^2),$$
the second equality obtained by expanding $(x_H+\delta/r_1)(y_H+\delta/r_2)$
to first order in $\delta$ and simplifying $x_H/r_2+y_H/r_1=-e_H^2 s$ (an
identity verifiable directly from the definitions of $x_H,y_H,e_H,s$).

**Step 2 -- both sides of $(\ast)$, to order $\varphi^4$.** With
$\delta=\varphi^2 q$ (so $\delta=O(\varphi^2)$, $\delta^2=O(\varphi^4)$) and
the standard expansions $\sin^2\varphi=\varphi^2-\varphi^4/3+O(\varphi^6)$,
$\cos\varphi=1-\varphi^2/2+\varphi^4/24+O(\varphi^6)$, writing
$e^2=e_H^2+\varphi^2 G+O(\varphi^4)$ (form to be confirmed), identity
$(\ast)$ -- $e^2\sin^2\Delta\theta=x^2+y^2-2xy\cos\Delta\theta$, with
$\cos\Delta\theta=-\cos\varphi$, $\sin\Delta\theta=\sin\varphi$ -- expands
term by term (direct algebraic computation, no numerical approximation) into
an equality of polynomials in $\varphi$ up to order $\varphi^4$. The
order-$\varphi^2$ term matches identically (no constraint on $q$, confirming
that $e_H^2$ is indeed the correct leading term); the order-$\varphi^4$ term
gives, after simplification,

$$\boxed{\,G(q) = s^2 q^2 + e_H^2\, s\, q + \frac{e_H^2}{4}\,}$$

*(The detail of the expansion -- a dozen lines of elementary algebra,
carrying forward Steps 1-2 -- is given in full in Appendix A; the identity
$d:=1/r_1-1/r_2=e_H\,s$, used for the cancellation at order $\varphi^2$, is
verified in one line from the definitions.)*

**Remark (positivity of the sub-question on $e^2$).** $G$ is a strictly
convex quadratic form in $q$ (leading coefficient $s^2>0$), with minimum
$G_{\min}=G(q^\ast_e)$ at $q^\ast_e=-e_H^2/(2s)$, equal to
$$G_{\min} = \frac{e_H^2}{4}-\frac{e_H^4}{4} = \frac{e_H^2(1-e_H^2)}{4} \;>\;0,$$
strictly positive since $0<|e_H|<1$ for $r_1\neq r_2$. This is a standalone
result (the minimum of $e^2$ along any path $p(\varphi)=p_H+\varphi^2q$
strictly increases away from $\varphi=0$) -- but **this is not yet** Lemma 3:
minimizing $e^2$ over $q$ at fixed $\varphi$ is not the same operation as
minimizing $\Delta v_T$ over $q$ at fixed $\varphi$ (the relationship between
the two, via Lemma 1, only applies at **fixed** $p$, whereas here $p$ varies
with $q$). The computation must be redone directly on $\Delta v_T$.

## 6. Expansion of $\Delta v_T$ near the Hohmann point

**A preliminary fact: the two passage points of the Hohmann transfer are its
apsides.** At $\varphi=0$ ($\delta=0$), $x_H=e_H$ and $y_H=-e_H$. Since
$e\cos\nu_k=\xi_k$ (Definition 1) and, at $\varphi=0$, $e=e_H$ (by continuity
of $e^2\to e_H^2$, Sec. 5), we get $\cos\nu_1=x_H/e_H=1$ and
$\cos\nu_2=y_H/e_H=-1$, i.e. $\nu_1=0$, $\nu_2=\pi$: the passage points are
exactly the perigee and apogee of the Hohmann transfer. In particular
$\sin\nu_1=\sin\nu_2=0$, so, by Definition 2, **the radial velocity component
vanishes there**: $v_{\theta,k,0}=v_{k,0}$ (subscript $0$ = value at
$\varphi=0$), an elementary and classical fact, here re-derived directly from
the framework set up in Definitions 1-2, not assumed without proof.

**Expansion of $v_k^2$.** By the vis-viva identity rewritten
($v_k^2=(\mu/p)e^2+2\mu/r_k-\mu/p$, Sec. 1) and the expansions
$1/p = (1/p_H)(1-\delta/p_H+O(\delta^2))$,
$e^2=e_H^2+\varphi^2G(q)+O(\varphi^4)$:

$$v_k^2 = v_{k,0}^2 + \varphi^2\, W(q) + O(\varphi^4), \qquad W(q):=\frac{\mu}{p_H}\left[G(q)-\frac{q\,(e_H^2-1)}{p_H}\right],$$

where $v_{k,0}^2 = \mu(2/r_k-1/a_H)$ (vis-viva on the Hohmann ellipse
itself, $a_H=(r_1+r_2)/2$). **Notable fact: $W(q)$ does not depend on $k$**
-- because $1/a=(1-e^2)/p$ is a single number for the whole ellipse, shared
by both passage points.

**Expansion of $\Delta v_k$.** $v_{\theta,k}=\sqrt{\mu p}/r_k =
v_{\theta,k,0}\left(1+\dfrac{\varphi^2 q}{2p_H}+O(\varphi^4)\right)$.
Substituting into $\Delta v_k^2=v_k^2+v_{c,k}^2-2v_{c,k}v_{\theta,k}$
(Definition 2) and using $v_{\theta,k,0}=v_{k,0}$:

$$\Delta v_k^2 = \Delta v_{k,0}^2 + \varphi^2\left[W(q) - \frac{v_{c,k}v_{k,0}}{p_H}\,q\right] + O(\varphi^4), \qquad \Delta v_{k,0}:=|v_{k,0}-v_{c,k}|,$$

where $\Delta v_{k,0}$ is the classical Hohmann impulse at the $k$-th point
($\Delta v_{1,0}+\Delta v_{2,0}=\Delta v_{\text{Hohmann}}$). Since
$\Delta v_{k,0}>0$ ($r_1\neq r_2$), take the square root by the same argument
as in Lemma 1:

$$\Delta v_k = \Delta v_{k,0} + \frac{\varphi^2}{2\Delta v_{k,0}}\left[W(q)-\frac{v_{c,k}v_{k,0}}{p_H}q\right] + O(\varphi^4).$$

**Summation.** Writing
$P:=\dfrac12\left(\dfrac1{\Delta v_{1,0}}+\dfrac1{\Delta v_{2,0}}\right)>0$ and
$Q:=\dfrac1{2p_H}\left(\dfrac{v_{c,1}v_{1,0}}{\Delta v_{1,0}}+\dfrac{v_{c,2}v_{2,0}}{\Delta v_{2,0}}\right)>0$:

$$\Delta v_T(p,\Delta\theta) \;=\; \Delta v_{\text{Hohmann}} \;+\; \varphi^2\,F(q) \;+\; O(\varphi^4), \qquad F(q):=P\,W(q)-Q\,q.$$

Substituting $W(q)$ and $G(q)$ (Sec. 5), $F$, like $G$, is a quadratic form
in $q$:

$$F(q) = a\,q^2+b\,q+c, \qquad a=\frac{\mu P s^2}{p_H},\quad c = \frac{\mu P e_H^2}{4 p_H},\quad b = \frac{\mu P e_H^2 s}{p_H}-\frac{\mu P(e_H^2-1)}{p_H^2}-Q.$$

$a>0$ (product of strictly positive quantities), so $F$ is strictly convex:
it has a unique minimum, attained at $q^\ast=-b/(2a)$, and

$$K(r_1,r_2) \;=\; \min_q F(q) \;=\; c-\frac{b^2}{4a}.$$

## 7. Closing Lemma 3

**What is established, entirely by hand.** The singularity at
$\Delta\theta=\pi$ is resolved: along the ansatz $p=p_H+\varphi^2q$,
$\Delta v_T$ admits a regular expansion in $\varphi$ whose order-$\varphi^2$
coefficient is a strictly convex, explicit quadratic form in $q$, with
coefficients $a,b,c$ given by closed formulas in terms of the classical
quantities of the Hohmann transfer alone ($v_{c,k}$, $v_{k,0}$,
$\Delta v_{k,0}$, $e_H$, $p_H$, $s$) -- not a shapeless expression of several
thousand terms. $K(r_1,r_2)=c-b^2/(4a)$ follows directly (completing the
square, an elementary computation).

**What remains -- reduced to a single concrete algebraic statement.**
Lemma 3 (and hence the result of Sec. 0) is now equivalent to the single
inequality
$$b^2 \;<\; 4ac,$$
with $a,b,c$ as above -- a strictly simpler problem to handle by hand (an
inequality between three explicit quantities, built from a limited number of
square roots: $\sqrt{r_1}$, $\sqrt{r_2}$, $\sqrt{a_H}$) than a direct
comparison of $\Delta v_T$ at two arbitrary angles. This is the next step.

*Confirmatory test (not constitutive of the proof).* The formula
$K=c-b^2/(4a)$ above was evaluated, for verification, at six radius ratios:
$R=1.2,\,1.5,\,2,\,3,\,5,\,10$. The six values obtained match, to the
significant digit, a table already established independently. This agreement
confirms that the algebraic reduction above (Sec. 5-6) is correct; it was
not used at any step to produce it.

## 8. The inequality $b^2<4ac$: two reductions, then the proof

**Lemma 4 (symmetry).** $K(r_1,r_2)=K(r_2,r_1)$.

*Proof.* For fixed $p,\Delta\theta$, the formulas of Definition 2 for
$\Delta v_1$ and $\Delta v_2$ depend, respectively, only on $(p,r_1)$ and
$(p,r_2)$ -- neither depends on the order in which the two points are
visited. Swapping the roles of $r_1$ and $r_2$ (traversing the same conic, of
the same $p$, in the opposite direction) therefore simply swaps the labels
"point 1" and "point 2", without changing either $\Delta v_1+\Delta v_2$ or
the angle $\Delta\theta$ separating them. Thus
$\Delta v_T(p,\Delta\theta;r_1,r_2)=\Delta v_T(p,\Delta\theta;r_2,r_1)$ as
functions of $p$, hence the minima over $p$ are equal, so
$K(r_1,r_2)=K(r_2,r_1)$. $\blacksquare$

**Lemma 5 (scale invariance of the sign).** For every $c>0$, $K(cr_1,cr_2)$
and $K(r_1,r_2)$ have the same sign.

*Proof.* Under $r_1\mapsto cr_1,\ r_2\mapsto cr_2$ at fixed $\mu$, the
vis-viva identity $v^2=\mu(2/r-1/a)$ shows that every characteristic velocity
of the problem ($v_{c,k}$, $v_{k,0}$, etc.) is multiplied by $1/\sqrt c$
(each term $2/r$ and $1/a$ being multiplied by $1/c$). The cost $\Delta v_T$,
a sum of such velocities, is therefore multiplied by $1/\sqrt c$ at fixed
$\Delta\theta$, for the same value of $p/c$ corresponding to a rescaled $p$
-- the minimization over $p$ commutes with this rescaling. $K(r_1,r_2)$, the
coefficient of $\varphi^2$ in this same $\Delta v_T$ (a quantity of the same
physical nature, a velocity), is therefore also multiplied by
$1/\sqrt c>0$: $K(cr_1,cr_2)=K(r_1,r_2)/\sqrt c$. Multiplying by a strictly
positive constant does not change the sign. $\blacksquare$

**Consequence of Lemmas 4-5.** To establish $K(r_1,r_2)>0$ for every
$r_1\neq r_2>0$, it suffices to do so for $r_1=1$, $\mu=1$, $r_2=R>1$ (the
case $R<1$ follows from Lemma 4, the case of an arbitrary $r_1$ from
Lemma 5).

**Rationalizing substitution.** For $R>1$, set $t:=\sqrt R>1$ and
$u:=\sqrt{1+t^2}=\sqrt{1+R}$. The formulas of Sec. 6 (with
$r_1=1,r_2=t^2,\mu=1$) become, after direct simplification,

$$v_{1,0}=\frac{t\sqrt2}{u},\quad v_{2,0}=\frac{\sqrt2}{tu},\quad v_{c,1}=1,\quad v_{c,2}=\frac1t,\quad p_H=\frac{2t^2}{u^2},\quad e_H=\frac{t^2-1}{u^2},\quad s=\frac{u^2}{t^2}.$$

For $t>1$ one has $v_{1,0}>v_{c,1}$ and $v_{c,2}>v_{2,0}$ (direct check),
hence the Hohmann impulses, without absolute value:

$$\Delta v_{1,0}=\frac{t\sqrt2-u}{u}, \qquad \Delta v_{2,0}=\frac{u-\sqrt2}{tu}.$$

**The identity that rationalizes everything.** One checks directly, by
difference of squares,

$$(t\sqrt2-u)(t\sqrt2+u) = 2t^2-u^2 = t^2-1, \qquad (u-\sqrt2)(u+\sqrt2)=u^2-2=t^2-1$$

(the second equality on each line using $u^2=1+t^2$). Since
$t>1 \Rightarrow t^2-1>0$, one can invert $\Delta v_{1,0}$ and $\Delta v_{2,0}$
by multiplying by the conjugate, never dividing by an expression of
uncertain sign:

$$\frac1{\Delta v_{1,0}} = \frac{u(t\sqrt2+u)}{t^2-1}, \qquad \frac1{\Delta v_{2,0}} = \frac{tu(u+\sqrt2)}{t^2-1}.$$

**Explicit computation of $a$, $b$, $c$.** Writing, for brevity,
$N:=2\sqrt2\,t+(1+t)u$. Substituting $1/\Delta v_{1,0}$, $1/\Delta v_{2,0}$
(above), $s=u^2/t^2$, $e_H=(t^2-1)/u^2$, $p_H=2t^2/u^2$ into $P$ (Sec. 6), one
obtains, after adding the two fractions and factoring out $u$:

$$P = \frac{u\,N}{2(t^2-1)}.$$

Substituted into $a=Ps^2/p_H$ and $c=Pe_H^2/(4p_H)$ (Sec. 6-7) -- two direct
substitutions, each requiring only multiplying fractions and simplifying
powers of $t$ and $u$:

$$a = \frac{u^7 N}{4t^6(t^2-1)}, \qquad\qquad c = \frac{(t^2-1)\,N}{16\,t^2\,u}.$$

Computing $Q$ (Sec. 6, substituting the rationalized $v_{c,1}v_{1,0}/\Delta v_{1,0}$
and $v_{c,2}v_{2,0}/\Delta v_{2,0}$, reducing to the common denominator
$t(t^2-1)$, then simplifying $u(t^2+1)=u\cdot u^2=u^3$) gives
$$Q = \frac{2u^2(t^3+1)+\sqrt2\,u^5}{4t^3(t^2-1)}.$$

Computing $b=Pe_H^2s/p_H - P(e_H^2-1)/p_H^2 - Q$ takes a bit more work. The
first two terms, both proportional to $uN$ (direct computation, as for $a$
and $c$), add to $uN(t^4+1)/\big(4t^4(t^2-1)\big)$ -- using here
$(t^2-1)^2+2t^2=t^4+1$. Subtracting $Q$ and reducing to the common
denominator $4t^4(t^2-1)$, the numerator is
$$uN(t^4+1) - 2tu^2(t^3+1) - \sqrt2\,t\,u^5.$$
Expanding $uN=2\sqrt2\,tu+(1+t)u^2$ and grouping the $u$ terms (coefficient
$2\sqrt2\,t(t^4+1)-\sqrt2\,tu^4$, where $u^4=(1+t^2)^2$ combines into
$\sqrt2\,t(t^2-1)^2$ after simplifying $2(t^4+1)-u^4=(t^2-1)^2$) and the $u^2$
terms (coefficient $(1+t)(t^4+1)-2t(t^3+1)=(t-1)^2(t+1)(t^2+1)$, a direct
factorization of the polynomial in $t$), the numerator becomes
$$\sqrt2\,t\,u\,(t^2-1)^2 + u^2(t-1)^2(t+1)(t^2+1).$$
Substituting $t^2+1=u^2$ in the second term (so $u^2(t^2+1)=u^4$) and
factoring out $(t-1)^2(t+1)$ from both terms (using
$(t^2-1)^2=(t-1)^2(t+1)^2$ for the first), the numerator becomes
$(t-1)^2(t+1)\,u\,\big[\sqrt2\,t(t+1)+u^3\big]$. Dividing by
$4t^4(t^2-1)=4t^4(t-1)(t+1)$, a factor $(t-1)(t+1)$ cancels, leaving:

$$\boxed{\,b = \frac{(t-1)\,u\,M}{4t^4}\,}, \qquad M:=\sqrt2\,t(t+1)+u^3.$$

**$4ac-b^2$ by difference of squares.** With $a=u^7N/\big(4t^6(t^2-1)\big)$
and $c=(t^2-1)N/(16t^2u)$, the factor $(t^2-1)$ cancels exactly between $a$
and $c$:

$$4ac = 4\cdot\frac{u^7N}{4t^6(t^2-1)}\cdot\frac{(t^2-1)N}{16t^2u} = \frac{u^6N^2}{16t^8}.$$

And $b^2 = (t-1)^2u^2M^2/(16t^8)$. So

$$4ac-b^2 = \frac{u^2}{16t^8}\Big[u^4N^2-(t-1)^2M^2\Big] = \frac{u^2}{16t^8}\big[u^2N-(t-1)M\big]\big[u^2N+(t-1)M\big],$$

a difference of two squares. Expanding each of the two factors (direct
computation from $N=2\sqrt2t+(1+t)u$, $M=\sqrt2t(t+1)+u^3$, and $u^2=1+t^2$ to
simplify $2u^2\mp(t^2-1)$ into $t^2+3$, respectively $3t^2+1$):

$$u^2N-(t-1)M = \sqrt2\,t(t^2+3)+2u^3, \qquad\qquad u^2N+(t-1)M = t\big[\sqrt2(3t^2+1)+2u^3\big].$$

Substituting back:

$$\boxed{\;4ac-b^2 \;=\; \frac{u^2\cdot t\,\big[\sqrt2\,t(t^2+3)+2u^3\big]\,\big[\sqrt2(3t^2+1)+2u^3\big]}{16\,t^8}\;=\;\frac{u^2\big[\sqrt2\,t(t^2+3)+2u^3\big]\big[\sqrt2(3t^2+1)+2u^3\big]}{16\,t^7}\;}$$

**Manifest positivity.** For $t>0$ (so $u=\sqrt{1+t^2}>0$), each of the three
factors in the numerator -- $u^2$, $\sqrt2\,t(t^2+3)+2u^3$,
$\sqrt2(3t^2+1)+2u^3$ -- is a sum of strictly positive terms, hence strictly
positive; the denominator $16t^7$ is too. The right-hand side is a quotient
of manifestly positive quantities:

$$4ac-b^2 \;>\;0 \quad\text{for every } t>0, \qquad\text{i.e.}\qquad b^2<4ac \quad\text{for every } R>0,\ R\neq1.$$

*(Confirmatory check: the boxed identity above was cross-checked by exact
symbolic computation against $4ac-b^2$ obtained independently from the
definitions of $a,b,c$ in Sec. 6-7 -- exact agreement, with no
approximation, over the entire symbolic computation, not just at a few
points. This test confirms the factorization above; it did not produce it --
every factorization step above (difference of squares, substitution
$u^2=1+t^2$, grouping of the coefficients of $t^4+1$ and
$(t-1)^2(t+1)(t^2+1)$) is reproducible by hand, term by term.)*

## 9. Theorem

> **Theorem (strict local minimum of the Hohmann transfer).** For all
> $r_1,r_2>0$ with $r_1\neq r_2$, there exists $\delta>0$ such that, for every
> $\Delta\theta\in(\pi-\delta,\pi+\delta)\setminus\{\pi\}$ (restricted to
> $(0,\pi]$ for the physical part of the domain),
> $$\min_{p} \Delta v_T(p,\Delta\theta) \;>\; \Delta v_{\text{Hohmann}}.$$
> More precisely,
> $\min_p \Delta v_T(p,\Delta\theta) = \Delta v_{\text{Hohmann}} + K(r_1,r_2)\,(\Delta\theta-\pi)^2+O\big((\Delta\theta-\pi)^4\big)$
> with $K(r_1,r_2)=c-b^2/(4a)>0$.

**Proof.** By Sec. 4bis, paths with $p-p_H=O(\varphi)$ are not competitive:
the minimum of $\Delta v_T(p,\pi-\varphi)$ over $p$, for $\varphi$ small
enough, is necessarily attained on the family
$p=p_H+\varphi^2q+O(\varphi^4)$. On this family, Sec. 5-6 give
$\Delta v_T = \Delta v_{\text{Hohmann}}+\varphi^2F(q)+O(\varphi^4)$ with $F$
strictly convex quadratic (Sec. 6-7); the coefficient of
$\varphi^2=(\pi-\Delta\theta)^2$ in $\min_q\Delta v_T$ is therefore
$K(r_1,r_2)=\min_qF(q)=c-b^2/(4a)$ (completing the square). By Lemmas 4-5, it
suffices to establish $K(1,R)>0$ for $R>1$; this is exactly the content of
Sec. 8. $\blacksquare$

*(Remark: this proof does not invoke Lemma 2 -- it proceeds directly through
the expansion of $\Delta v_T$ itself on the competitive family, not through
the sign of $\partial(e^2)/\partial\Delta\theta$ at fixed $p$. Lemma 2
remains a correct and useful result in its own right -- for monotonicity at
fixed $p$, or as a building block of a possible future envelope argument
covering all of $(0,\pi)$ -- but it is not a necessary link in the chain
leading to this theorem.)*

**What this theorem establishes, and does not claim to establish.** This is
a **local** result -- a neighborhood $(\pi-\delta,\pi+\delta)$, not the whole
interval $(0,\pi)$. Prussing's theorem (1992) establishes that
$(p_{\text{Hohmann}},\pi)$ is the **global** minimum in the $(p,e)$ plane,
without ever fixing $\Delta\theta$; the present theorem establishes the
strict curvature, in $\Delta\theta$, of the cost already minimized over $p$
-- a different question, which Prussing's proof does not cover and could not
cover by its method. The two results complement each other: neither replaces
the other. Together, they characterize the Hohmann point as both a global
minimum *and* a strict local minimum in angle -- a property of shape (the
"valley" has no false plateau around its lowest point), not merely of
position. Monotonicity over the whole interval $(0,\pi)$ -- beyond this local
neighborhood -- remains an open, distinct question, not addressed here.
Likewise, this result concerns exclusively two-impulse transfers: beyond
$R=r_2/r_1\approx11.94$, a three-impulse bi-elliptic transfer can cost less
(Prussing, 1992) -- outside the scope of this document, which by construction
considers only a single Lambert arc (Definition 1).

---

## Summary of the proof chain

Lemma 1 (Sec. 2) and Lemma 2 (Sec. 3, an immediate corollary, not used in the
final proof of the theorem -- see the remark in Sec. 9 -- but a useful result
in its own right) set up the monotonicity of the cost in $e^2$ at fixed $p$.
Sec. 4bis establishes that the family $p-p_H=O(\varphi^2)$ is the only
competitive one -- not an ansatz, a consequence of Lemma 1 and parity.
"Lemma 3" (characterization of $p^*$ near $\pi$) is then resolved by hand
(Sec. 5-8), via Lemma 4 (symmetry) and Lemma 5 (scale invariance of the
sign), which reduce the problem to $r_1=1,R>1$, followed by an explicit
rationalization ($t=\sqrt R$, $u=\sqrt{1+t^2}$) and a difference-of-squares
factorization giving $4ac-b^2$ as a **product of three manifestly positive
factors**, divided by $16t^7>0$ -- positivity visible by inspection, with no
sum of polynomials or numerical evaluation needed. Every step of this
factorization (computation of $P,Q,a,b,c$, grouping into
$(t-1)^2(t+1)(t^2+1)$, substitution $t^2+1=u^2$, final difference of
squares) is written out in the main text, reproducible term by term by hand.
The central theorem (Sec. 9) then follows: the Hohmann transfer is a strict
local minimum, for every $r_1\neq r_2>0$, by an entirely hand-written proof.

## References

Prussing, J. E., "Simple Proof of the Global Optimality of the Hohmann
Transfer," *Journal of Guidance, Control, and Dynamics*, vol. 15, no. 4
(1992), pp. 1037-1038.

---

## Appendix A. Detail of the Sec. 5 expansion (identity $(\ast)$ to order $\varphi^4$)

Write, for brevity, $d:=\dfrac1{r_1}-\dfrac1{r_2}$ (not to be confused with
$s=1/r_1+1/r_2$). First, from the definitions $x_H=e_H$, $y_H=-e_H$,
$e_H=(r_2-r_1)/(r_1+r_2)$ and $s=(r_1+r_2)/(r_1r_2)$:

$$d = \frac{r_2-r_1}{r_1 r_2} = \frac{r_2-r_1}{r_1+r_2}\cdot\frac{r_1+r_2}{r_1r_2} = e_H\, s. \tag{A.1}$$

**Left-hand side.** With $\sin\varphi=\varphi-\varphi^3/6+O(\varphi^5)$, so
$\sin^2\varphi=\varphi^2-\varphi^4/3+O(\varphi^6)$, and
$e^2=e_H^2+\varphi^2G+O(\varphi^4)$ (form to be confirmed):

$$e^2\sin^2\varphi = (e_H^2+\varphi^2G)\Big(\varphi^2-\frac{\varphi^4}{3}\Big)+O(\varphi^6) = e_H^2\varphi^2 + \Big(G-\frac{e_H^2}{3}\Big)\varphi^4+O(\varphi^6). \tag{A.2}$$

**Right-hand side -- term $x^2+y^2$.** With $x=x_H+\delta/r_1$,
$y=y_H+\delta/r_2$, $\delta=\varphi^2q$:

$$x^2+y^2 = \underbrace{x_H^2+y_H^2}_{=2e_H^2} + 2\delta\Big(\frac{x_H}{r_1}+\frac{y_H}{r_2}\Big) + \delta^2\Big(\frac1{r_1^2}+\frac1{r_2^2}\Big).$$

The coefficient of $\delta$ is $x_H/r_1+y_H/r_2 = e_H/r_1-e_H/r_2 = e_H d$ (by
(A.1), $=e_H^2 s$). With $\delta=\varphi^2q$ (so $\delta^2=\varphi^4q^2$):

$$x^2+y^2 = 2e_H^2 + 2\varphi^2 q\, e_H^2 s + \varphi^4 q^2\Big(\frac1{r_1^2}+\frac1{r_2^2}\Big)+O(\varphi^6). \tag{A.3}$$

**Right-hand side -- term $xy$.**
$xy = x_Hy_H+\delta\big(\frac{x_H}{r_2}+\frac{y_H}{r_1}\big)+\delta^2/(r_1r_2)$,
with $x_Hy_H=-e_H^2$ and, by a computation analogous to (A.1),
$$\frac{x_H}{r_2}+\frac{y_H}{r_1} = e_H\Big(\frac1{r_2}-\frac1{r_1}\Big) = -e_H d = -e_H^2 s.$$
Hence
$$xy = -e_H^2 - \varphi^2 q\, e_H^2 s + \varphi^4\frac{q^2}{r_1r_2}+O(\varphi^6). \tag{A.4}$$

**Right-hand side -- term $-2xy\cos\Delta\theta$.** With
$\cos\Delta\theta=-\cos\varphi=-1+\varphi^2/2-\varphi^4/24+O(\varphi^6)$ and
(A.4):

$$-2xy\cos\Delta\theta = 2xy - xy\,\varphi^2 + \frac{xy}{12}\varphi^4+O(\varphi^6)$$
$$= \Big[-2e_H^2-2\varphi^2qe_H^2s+2\varphi^4\frac{q^2}{r_1r_2}\Big] + \Big[e_H^2\varphi^2+\varphi^4qe_H^2s\Big] + \Big[-\frac{e_H^2}{12}\varphi^4\Big]+O(\varphi^6)$$
$$= -2e_H^2 + \varphi^2\big(e_H^2-2qe_H^2s\big) + \varphi^4\Big(\frac{2q^2}{r_1r_2}+qe_H^2s-\frac{e_H^2}{12}\Big)+O(\varphi^6). \tag{A.5}$$

**Sum (A.3)+(A.5).** The constant term cancels ($2e_H^2-2e_H^2=0$,
consistent with $e^2\to e_H^2$). Grouping:

$$\text{RHS} = \varphi^2\Big[e_H^2+2qe_H^2s-2qe_H^2s\Big] + \varphi^4\Big[q^2\Big(\frac1{r_1^2}+\frac1{r_2^2}+\frac2{r_1r_2}\Big)+qe_H^2s-\frac{e_H^2}{12}\Big]+O(\varphi^6)$$

The $\varphi^2$ bracket simplifies immediately to $e_H^2$ (the two $q$ terms
cancel -- confirming that $e_H^2$ is indeed the leading term, with no
constraint on $q$). The $\varphi^4$ bracket uses
$\frac1{r_1^2}+\frac1{r_2^2}+\frac2{r_1r_2}=\big(\frac1{r_1}+\frac1{r_2}\big)^2=s^2$:

$$\text{RHS} = e_H^2\varphi^2 + \Big[q^2s^2+qe_H^2s-\frac{e_H^2}{12}\Big]\varphi^4+O(\varphi^6). \tag{A.6}$$

**Identification.** Equating (A.2) and (A.6), the order-$\varphi^2$ term
coincides identically; the order-$\varphi^4$ term gives
$$G-\frac{e_H^2}{3} = q^2s^2+qe_H^2s-\frac{e_H^2}{12} \quad\Longrightarrow\quad G(q) = q^2s^2+qe_H^2s+\frac{e_H^2}{4}$$
(using $-\tfrac1{12}+\tfrac13=\tfrac14$), which is exactly the boxed formula
in Sec. 5. $\blacksquare$

## Appendix B. Detail of the Sec. 8 computation ($P$, $Q$, $a$, $b$, $c$ in terms of $t$)

With $t=\sqrt R>1$, $u=\sqrt{1+t^2}$, $\Delta v_{1,0}=(t\sqrt2-u)/u$,
$\Delta v_{2,0}=(u-\sqrt2)/(tu)$, rationalizing by the conjugate (Sec. 8)
gives
$$\frac1{\Delta v_{1,0}}=\frac{u(t\sqrt2+u)}{t^2-1}, \qquad \frac1{\Delta v_{2,0}}=\frac{tu(u+\sqrt2)}{t^2-1}.$$

**Computation of $P$.**
$2P=\dfrac1{\Delta v_{1,0}}+\dfrac1{\Delta v_{2,0}} = \dfrac{u(t\sqrt2+u)+tu(u+\sqrt2)}{t^2-1}$.
The numerator expands:
$u(t\sqrt2+u)+tu(u+\sqrt2)=t\sqrt2\,u+u^2+tu^2+t\sqrt2\,u = 2t\sqrt2\,u+(1+t)u^2$.
Hence

$$P = \frac{u\big[2t\sqrt2+(1+t)u\big]}{2(t^2-1)}. \tag{B.1}$$

**Computation of $Q$.** From
$v_{c,1}v_{1,0}/\Delta v_{1,0}=t\sqrt2(t\sqrt2+u)/(t^2-1)$ and
$v_{c,2}v_{2,0}/\Delta v_{2,0}=\sqrt2(u+\sqrt2)/\big(t(t^2-1)\big)$ (direct
computation from the rationalized expressions above and from
$v_{c,1}=1,v_{c,2}=1/t,v_{1,0}=t\sqrt2/u,v_{2,0}=\sqrt2/(tu)$), reducing to
the common denominator and using $u^2=1+t^2$ to simplify $u(t^2+1)=u^3$:

$$v_{c,1}v_{1,0}/\Delta v_{1,0}+v_{c,2}v_{2,0}/\Delta v_{2,0} = \frac{\sqrt2\big[\sqrt2(t^3+1)+u^3\big]}{t(t^2-1)},$$

then, with $p_H=2t^2/u^2$,

$$Q = \frac{2u^2(t^3+1)+\sqrt2\,u^5}{4t^3(t^2-1)}. \tag{B.2}$$

**From $P,Q$ to $a,b,c$, then to $4ac-b^2$.** The rest of the computation --
substituting (B.1)-(B.2) into $a=Ps^2/p_H$, $c=Pe_H^2/(4p_H)$,
$b=Pe_H^2s/p_H-P(e_H^2-1)/p_H^2-Q$, then factoring $4ac-b^2$ by difference of
squares -- is written out, step by step, in the main text (Sec. 8, starting
at "Explicit computation of $a,b,c$"), not relegated to this appendix: each
simplification (factoring out $(t-1)^2(t+1)(t^2+1)$, substituting
$t^2+1=u^2$, the difference of squares $u^2N\mp(t-1)M$) is shown there, not
merely asserted.
