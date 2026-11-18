# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
"""
Exact algebraic certificate for Phi > 0.

Variables:
    0 < w < 1, z > 0, u > 0, where
        w = R^{-1/2}, z = sqrt(p), u = q/z,
        q = cot(phi/2).

The proof uses:
  H = S1 + S2
    = d1*u^2 + d0 + dminus/u^2

  K/2 = A_lambda + 3*q^2/z^4
      = c1*u^2 + c0 + cminus/u^2

and Cauchy-Schwarz:
  (K/2) H >=
    (sqrt(c1*dminus) + sqrt(c0*d0) + sqrt(cminus*d1))^2.

The three square-root terms are then bounded from below.
No numerical argument is used.
"""

import sympy as sp

w, z, u = sp.symbols("w z u", positive=True)

a = 1 - w**2
b = 1 + w**2
A0 = 1 + w**4
m = 1 - w**3

q = u*z

lam_minus = (b*q - a/q)/2
A_lambda = 1 + lam_minus**2

S1 = sp.expand((z - 1)**2 + (lam_minus*z - q/z)**2)
lam_plus = (b*q + a/q)/2
S2 = sp.expand(
    w**2*(w*z - 1)**2 + (lam_plus*z - q/z)**2
)

H = sp.factor(S1 + S2)
K_over_2 = sp.factor(A_lambda + 3*q**2/z**4)

R0 = (z - 1)**2 + w**2*(w*z - 1)**2

d1 = sp.factor((b*z**2 - 2)**2/2)
d0 = sp.factor(R0)
dm = sp.factor(a**2/2)

c1 = sp.factor(b**2*z**2/4 + 3/z**2)
c0 = sp.factor(A0/2)
cm = sp.factor(a**2/(4*z**2))

assert sp.factor(H - (d1*u**2 + d0 + dm/u**2)) == 0
assert sp.factor(K_over_2 - (c1*u**2 + c0 + cm/u**2)) == 0

print("Exact decompositions verified.")

# Cauchy certificate:
# vectors:
#   X = (sqrt(c1)u, sqrt(c0), sqrt(cm)/u)
#   Y = (sqrt(dm)/u, sqrt(d0), sqrt(d1)u)
#
# Hence K/2 * H >= L^2, with
# L = sqrt(c1*dm)+sqrt(c0*d0)+sqrt(cm*d1).

T1 = sp.factor(
    a/(2*sp.sqrt(2)*z) * sp.sqrt(b**2*z**4 + 12)
)
T0 = sp.sqrt(A0*R0/2)
T3 = sp.factor(
    a*sp.Abs(b*z**2 - 2)/(2*sp.sqrt(2)*z)
)

print("Cauchy terms:")
print("T1 =", T1)
print("T0 =", T0)
print("T3 =", T3)

# Algebraic identity for the radial quadratic:
z0 = sp.factor((1 + w**3)/A0)
R_complete = sp.factor(
    A0*(z - z0)**2 + w**2*(1 - w)**2/A0
)
assert sp.factor(R0 - R_complete) == 0
print("Radial completion of the square verified.")

print("\nConsequences used in the proof:")
print("T0 >= w*(1-w)/sqrt(2).")
print("T1+T3 >= a*sqrt(b), via the lemma below.")

# Verify the one-variable lemma
X = sp.symbols("X", positive=True)

# Lemma:
# sqrt(X^2+12) + |X-2| >= 2*sqrt(2X)
#
# It suffices to square after moving |X-2| to the other side.
# The two branches reduce to the following exact factorizations
# after X=t^2:
t = sp.symbols("t", nonnegative=True)

low = sp.factor(
    t**4 + 12 - (t**2 - 2 + 2*sp.sqrt(2)*t)**2
)
high = sp.factor(
    t**4 + 12 - (2*sp.sqrt(2)*t - (t**2 - 2))**2
)

print("\nOne-variable lemma factorisations:")
print("0 < X <= 2:", low)
print("X >= 2:", high)

# The RHS of the squared inequality is nonnegative precisely
# in the nontrivial subranges; the displayed factorizations give
# the sign there because t^2-2 changes sign at sqrt(2).

# Final elementary comparison:
# m - a*sqrt(b)
# = (1-w) * [1+w+w^2 - (1+w)*sqrt(1+w^2)]
# = (1-w)*w^2 / denominator.
#
# The denominator > sqrt(2)*w, so the deficit is < (1-w)*w/sqrt(2).
#
# The symbolic rationalization identity:
den = 1 + w + w**2 + (1 + w)*sp.sqrt(b)
deficit = sp.factor(
    (1 + w + w**2) - (1 + w)*sp.sqrt(b)
)
rationalized = sp.simplify(
    deficit * den - w**2
)

print("\nRationalization check:")
print("(deficit)*den - w^2 =", rationalized)

assert rationalized == 0

print("\nCertificate summary:")
print("L >= a*sqrt(b) + w*(1-w)/sqrt(2) > m.")
print("Therefore (K/2)*(S1+S2) > m^2.")
print("At every stationary point, C_zz > 0.")
print("Since C_zz = 4 p C_pp at a stationary point, H2 follows.")
print("Moreover every stationary point is a strict local minimum, hence H1 follows.")
