# Theorem 4 — Proof Map and Global Certification Strategy

## 1. Objective

The purpose of this document is to describe how the local analytical arguments and validated numerical computations combine to establish the global sign pattern asserted in Theorem 4.

We consider the cost function C(p,φ) and its physical stationary branch Γ_phys, defined by the stationary equation

    C_p(p,φ) = 0

together with the physical constraints

    0 < p < 4/3,   0 < φ < π.

The existence, uniqueness, and interval enclosures of the relevant stationary solution are established in the certified regions described below; the constraints above alone do not distinguish Γ_phys from any other solution of C_p = 0 that might lie in this domain.

The quantity of interest is the determinant of the Hessian matrix of C:

    H(p,φ) = C_φφ C_pp − C_pφ².

We use

    s = √p,   s_H = 2/√3,

as the branch parameter in this proof. This notation is deliberately distinct from the parameter w = tan(φ/2) used elsewhere in the article.

The objective is to prove that there is exactly one point s_c on the physical stationary branch at which H vanishes, with

    H(s) < 0   for 0 < s < s_c,
    H(s) > 0   for s_c < s < s_H.

The proof is global: it includes both singular endpoint regimes and all junctions between the interval regions.

## 2. Mathematical formulation

The cost function is expressed through two positive radicands:

    C(p,φ) = √S1(p,φ) + √S2(p,φ).

For the R=2 model, introduce

    x = cos φ,   A0 = p − 1,   B0 = p/2 − 1,
    N = A0 x + B0,   E = A0² + N² / (1 − x²).

Then

    S1 = 3 + (E − 1)/p − 2√p,
    S2 = 3/2 + (E − 1)/p − √(p/2).

The positivity of S1 and S2 is checked on the interval domains used in the proof. This is necessary for the differentiability and validity of the square-root derivative formulas.

For q, r in {p, x}, the derivatives of C are computed from

    C_q = Σ_{i=1,2} (S_i)_q / (2√S_i),

and

    C_qr = Σ_{i=1,2} [ (S_i)_qr / (2√S_i) − (S_i)_q (S_i)_r / (4 S_i^{3/2}) ].

The angular derivatives are recovered using

    C_φφ = (1 − x²) C_xx − x C_x,
    C_pφ = −sin φ · C_px.

Consequently,

    H = [ (1 − x²) C_xx − x C_x ] C_pp − (1 − x²) C_px².

All interval evaluations must use a consistent convention for these transformations and their signs.

## 3. Why the proof is divided into regions

A direct interval mesh in the original variables is not well suited to the entire stationary branch.

Near the lower endpoint, the physical branch approaches a degenerate configuration with x → −1. Near the Hohmann endpoint, the angular variable approaches φ = 0, and the natural endpoint coordinate is singular in the original angular derivatives.

The proof therefore uses different regularizations near the two endpoints and a conventional interval mesh in the interior.

The six numbered directories correspond to the following six stages:

1. lower endpoint regularization;
2. lower bridge;
3. central mesh and critical-point isolation;
4. Hohmann endpoint regularization;
5. upper bridge;
6. global audit.

The proof is complete only when the local certificates and all their junctions have been checked.

## 4. Lower endpoint: s → 0

**Directory:** `01_lower_endpoint/`

The lower endpoint is regularized by introducing

    x = −1 + s² Y,   s = √p.

This substitution resolves the degeneracy in 1 − x² and provides a bounded auxiliary variable Y on the relevant endpoint box.

The stationary equation is regularized by defining

    G(s, Y) = s² C_p(s², −1 + s² Y).

The endpoint analysis establishes the existence and uniqueness of the physical stationary solution in the prescribed Y-interval. It also establishes the sign of the appropriately rescaled Hessian determinant.

The branch has the expansion

    x(s) = −1 + s²/2 − (1/3 + √2/6) s³ + (2 + √2)/4 s⁴ + O(s⁵).

Equivalently, since x = −1 + s² Y, the auxiliary variable satisfies

    Y(s) = 1/2 − (1/3 + √2/6) s + (2 + √2)/4 s² + O(s³).

The sign certificate uses the rescaled quantity

    P(s) = s⁶ H(s).

The endpoint calculation establishes

    (1/2) P''(0) = −(9 + 4√2)/32 < 0.

Together with the interval enclosure of P'' on the endpoint box, this gives

    H(s) < 0

throughout the certified lower endpoint region.

The corresponding computational certificate is produced by `certify_Ppp_v2.py`.

## 5. Lower bridge

**Directory:** `02_lower_bridge/`

The lower endpoint regularization is connected to the central interval mesh by a validated interval bridge.

The reported bridge covers

    0.00015 ≤ s ≤ 0.15

using 464 accepted interval tiles.

On every tile, the computation verifies the relevant stationary-branch inclusion, the required radicand positivity, and the strict negativity of H.

The junction at the lower endpoint is checked through the common auxiliary-variable enclosure and the endpoint signs of the regularized stationary equation.

The junction with the central mesh is checked through overlapping enclosures and interval inclusion for the stationary branch.

The pilot mesh, if retained, is a redundant check and is not part of the minimal global coverage argument.

The main certificate is associated with `full_mesh_bridge.py`.

## 6. Central interval mesh and critical point

**Directory:** `03_central_mesh/`

The central region is treated in the original variables, with interval enclosures of the stationary branch and the Hessian determinant.

The reported central computation divides the region into two parts separated by a small box containing the critical point.

### 6.1. Lower central part

The lower central mesh covers the interval from s = 0.15 to the lower boundary s_c⁻ of the critical-point box.

The computation uses 2052 accepted tiles and establishes

    H < 0

throughout this region.

### 6.2. Critical point isolation

The critical point is characterized by the system

    C_p(p,φ) = 0,   H(p,φ) = 0.

A Krawczyk operator is used to prove existence and uniqueness of a solution in the prescribed critical-point box.

The reported numerical approximation is

    s_c   ≈ 0.3069808139018763087165915174,
    p_c   ≈ 0.0942372201038584147638973774,
    φ_c   ≈ 2.8688160463235594826404957334  radians.

The rigorous result is the interval enclosure and the strict Krawczyk inclusion, not the decimal approximation itself.

The certified critical interval in the branch parameter s = √p is

    s_c⁻ = 0.306980651025204542554804388093645823472607040819467619479604,
    s_c⁺ = 0.306980976778461656460911810468050558140482032920939162941821.

These endpoints are the outward-rounded bounds used in the global coverage argument. The Krawczyk certificate establishes existence and uniqueness of the solution of

    (C_p(p,φ), H(p,φ)) = (0, 0)

in the certified box in the original (p,φ) variables. The interval above records its projection onto the branch parameter s = √p; it is not, by itself, the Krawczyk box. Both objects — the (p,φ) Krawczyk box and its s-projection — must be preserved in the certificate file.

### 6.3. Upper central part

The upper central mesh covers the interval from the upper boundary s_c⁺ of the critical-point box to the central upper boundary

    s = 1.144701.

The reported computation uses 3090 accepted tiles and establishes

    H > 0

throughout this region.

The two strict-sign meshes, together with the unique zero isolated by the Krawczyk computation, exclude any additional zero of H in the central region.

The central certificate is produced by `mesh_certificate_R2.py`.

## 7. Hohmann endpoint: p → 4/3

**Directory:** `04_hohmann_endpoint/`

At the upper endpoint,

    p_H = 4/3,   s_H = 2/√3,   φ → 0.

The endpoint is regularized by introducing

    t = sin²(φ/2),   δ = p − 4/3.

The stationary branch has the expansion

    p(t) = 4/3 + a1 t + O(t²),

where

    a1 = (88√3 − 56√6 − 24√2 − 28) / 621.

The endpoint radicands are strictly positive:

    S1_H = 7/3 − 4√3/3 > 0,
    S2_H = 5/6 − √6/3 > 0.

The endpoint analysis establishes the asymptotic relation

    t H(p(t), t) → L,

with

    L = (105 + 81√2 + 63√3 + 45√6) / 256 > 0.

Thus the limiting behavior is consistent with strict positivity of H near the Hohmann endpoint.

The proof also includes a finite-t interval certificate. The limiting formula alone is not used as a substitute for a uniform sign enclosure on a nonzero interval.

The endpoint computation is associated with `hohmann_certificate_v3.py`; the exact conversion of the endpoint bound to the s-parameter is checked by `hohmann_w_bound_exact.py`.

## 8. Upper bridge

**Directory:** `05_upper_bridge/`

The upper bridge connects the central mesh to the Hohmann endpoint region.

It is parameterized using t, and its interval bounds are mapped back to the s = √p parameter using the certified stationary-branch enclosure and the monotonicity and conversion checks.

The reported upper bridge contains 98 accepted tiles and establishes strict positivity of H throughout its domain.

Two junctions are essential:

1. the junction at t = 10⁻⁴ with the local Hohmann endpoint certificate;
2. the junction with the central mesh.

The first is checked by strict inclusion of the bridge auxiliary-variable enclosure in the local endpoint box, together with the relevant existence and uniqueness statements.

The second is checked through overlapping enclosures in (p,φ) and a strict interval sign for C_pφ on the common overlap. This establishes the required identification of the stationary-branch solutions at the junction.

The exact conversion of the lower s-bound of the upper bridge is also checked, including the strict inequality

    s_- < s_H.

The principal scripts are `hohmann_full_mesh.py`, `hohmann_raccord_t1e4.py`, and `hohmann_central_join.py`.

## 9. Global coverage

The reported global coverage is summarized below.

| Region             | Parameter interval      | Method                               | Certified sign |
| ------------------- | ------------------------- | --------------------------------------- | ----------------- |
| Lower endpoint     | 0 < s ≤ 0.00015         | Regularized endpoint certificate     | H < 0          |
| Lower bridge       | 0.00015 ≤ s ≤ 0.15      | Validated interval mesh              | H < 0          |
| Central lower part | 0.15 ≤ s ≤ s_c⁻         | Central interval mesh                | H < 0          |
| Critical box       | s_c⁻ ≤ s ≤ s_c⁺         | Krawczyk isolation of (C_p, H) = 0   | Unique zero    |
| Central upper part | s_c⁺ ≤ s ≤ 1.144701     | Central interval mesh                | H > 0          |
| Upper bridge       | 1.144701 ≤ s ≤ s_-      | Validated mesh in t, mapped to s     | H > 0          |
| Hohmann endpoint   | s_- ≤ s < s_H           | Regularized endpoint certificate     | H > 0          |

Here s_c⁻ and s_c⁺ denote the certified lower and upper bounds of the critical-point enclosure. The value s_- is the certified lower s-bound obtained from the upper bridge endpoint conversion.

The numerical value 1.144701 is the central/upper-bridge junction in the s = √p parameter. The upper bridge itself is computed in the t-parameter; the table records its certified image in s.

The upper bridge endpoint t_TARGET is a finite decimal endpoint used in the computation. It should be treated as the exact decimal value specified in the certificate, not as an exact symbolic representation of an unspecified limiting point.

The global audit must verify that these regions cover the full open physical interval

    0 < s < s_H,

with compatible branch enclosures at every junction and no uncovered subinterval.

## 10. Interpretation of the result

Since the stationary branch satisfies

    C_p = 0

and C_pp > 0 holds throughout the certified regions (C_pp, or the equivalent monotonicity quantity t·C_pp on the Hohmann side, is part of the existence-and-uniqueness argument in each region; its consolidation into a single global statement across all seven regions is part of the global audit, see Section 11), the second derivative of the cost restricted to the stationary branch is

    d²/dφ² C(p*(φ), φ) = (C_φφ C_pp − C_pφ²) / C_pp = H / C_pp.

Therefore, the sign change of H is also the sign change of the curvature of the cost restricted to the stationary branch.

In particular:

- for s < s_c, the Hessian determinant is negative and the Hessian is indefinite;
- at s = s_c, the Hessian determinant vanishes;
- for s > s_c, the Hessian determinant is positive; together with C_pp > 0, the Hessian is positive definite.

The sign pattern of H is the primary conclusion of the global certification. The interpretation in terms of Hessian definiteness additionally requires the global positivity statement C_pp > 0, which is to be consolidated by the audit described in Section 11; on the stationary branch, the curvature formula H / C_pp only requires C_pp ≠ 0, whereas concluding that the Hessian is positive definite where H > 0, or that the branch is a local minimum in p, requires the stronger C_pp > 0.

Once this global positivity check has passed, the unique zero of H marks the transition from an indefinite Hessian to a positive-definite Hessian along the physical stationary branch. Equivalently, it marks the sign change of the curvature of the cost restricted to that branch.

## 11. Audit and reproducibility

The global conclusion depends on the following checks:

- exact symbolic identities used in the endpoint regularizations;
- positivity of the radicands on every certified interval domain;
- consolidation of C_pp > 0 into a single global statement spanning all seven regions, rather than a fact used only piecewise inside each region's existence argument. The quantity actually certified, and the condition required of it, differ by region:

  | Region           | Certified quantity          | Required bound        |
  | ----------------- | ---------------------------- | ------------------------ |
  | Lower endpoint    | C_pp, or a rescaled form with a strictly positive rescaling factor | lower bound > 0 |
  | Lower bridge      | C_pp                         | lower bound > 0       |
  | Central lower     | C_pp                         | lower bound > 0       |
  | Critical box      | C_pp, over the full relevant box | lower bound > 0   |
  | Central upper     | C_pp                         | lower bound > 0       |
  | Upper bridge      | t·C_pp, with t > 0 on the region | lower bound > 0   |
  | Hohmann endpoint  | t·C_pp, for t > 0            | lower bound > 0       |

  A positive lower bound for t·C_pp establishes C_pp > 0 only where the certificate also establishes t > 0 on the domain considered; at the limiting point t = 0 this must be argued by a limit or by the regularized expansion, not by dividing by t directly;
- correct signs and factors in all variable transformations;
- existence and uniqueness of the stationary branch on the relevant boxes;
- strict sign enclosures for H on all noncritical regions;
- strict Krawczyk inclusion for the unique critical point;
- overlap and branch-identification checks at all regional junctions;
- exact conversion and ordering of the global parameter bounds.

The scripts in `06_global_audit/` are intended to collect these checks into a single global verification. In particular, the global C_pp audit must: (1) read the certified bounds, not the floating-point center values; (2) verify the rescaling factors and their sign; (3) verify the exact domains covered; (4) verify that the seven regions jointly cover the full interval 0 < s < s_H; (5) reject the global certificate as soon as one region fails to supply an exploitable, strictly positive lower bound. This audit must work directly from the archived certificates and their recorded bounds, without re-executing the producing scripts and without relying on their printed conclusions.

The archived certificates, their generating scripts, software versions, and SHA-256 checksums should be linked through `data/certificate_manifest.json`.

For details on running the computations and verifying the archived results, see `REPRODUCIBILITY.md`.

## 12. Scope of the certification

The certification described here concerns the R=2 model and the physical stationary branch specified above.

The numerical approximations are provided for readability. The mathematical conclusions rely on the exact identities, interval enclosures, strict inclusion tests, and global coverage checks archived with the computational proof.

The proof should be considered reproducible from a given repository release only when all certificate files, their generating scripts, the dependency environment, and the global audit are present in that release.

## Copyright

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
