# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
import sympy as sp
import time, pickle

w, Y, u1, u2 = sp.symbols('w Y u1 u2', positive=True)
with open('/tmp/P_built.pkl','rb') as f:
    P = pickle.load(f)

S1 = (8*Y**2*w**3 - 12*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 16*Y*w + 12*Y + 1) / (4*Y*(2 - Y*w**2))
S2 = (2*sp.sqrt(2)*Y**2*w**3 - 6*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 4*sp.sqrt(2)*Y*w + 1) / (4*Y*(2 - Y*w**2))

t0=time.time()
S1_w = sp.diff(S1,w); S1_Y = sp.diff(S1,Y)
S1_ww = sp.diff(S1,w,w); S1_wY = sp.diff(S1,w,Y); S1_YY = sp.diff(S1,Y,Y)
S2_w = sp.diff(S2,w); S2_Y = sp.diff(S2,Y)
S2_ww = sp.diff(S2,w,w); S2_wY = sp.diff(S2,w,Y); S2_YY = sp.diff(S2,Y,Y)
print(f"[{time.time()-t0:.1f}s] S1,S2 derivatives done")

# u1=sqrt(S1): standard formulas
u1_w = S1_w/(2*u1); u1_Y = S1_Y/(2*u1)
u1_ww = S1_ww/(2*u1) - S1_w**2/(4*u1**3)
u1_wY = S1_wY/(2*u1) - S1_w*S1_Y/(4*u1**3)
u1_YY = S1_YY/(2*u1) - S1_Y**2/(4*u1**3)
u2_w = S2_w/(2*u2); u2_Y = S2_Y/(2*u2)
u2_ww = S2_ww/(2*u2) - S2_w**2/(4*u2**3)
u2_wY = S2_wY/(2*u2) - S2_w*S2_Y/(4*u2**3)
u2_YY = S2_YY/(2*u2) - S2_Y**2/(4*u2**3)

t0=time.time()
P_w = sp.diff(P,w); P_Y = sp.diff(P,Y); P_u1 = sp.diff(P,u1); P_u2 = sp.diff(P,u2)
print(f"[{time.time()-t0:.1f}s] P first partials (4) done, ops={[sp.count_ops(e) for e in (P_w,P_Y,P_u1,P_u2)]}")

t0=time.time()
P_ww = sp.diff(P_w,w); P_wY = sp.diff(P_w,Y)
P_YY = sp.diff(P_Y,Y)
P_w_u1 = sp.diff(P_w,u1); P_w_u2 = sp.diff(P_w,u2)
P_Y_u1 = sp.diff(P_Y,u1); P_Y_u2 = sp.diff(P_Y,u2)
P_u1u1 = sp.diff(P_u1,u1); P_u1u2 = sp.diff(P_u1,u2); P_u2u2 = sp.diff(P_u2,u2)
print(f"[{time.time()-t0:.1f}s] P second partials (9) done")

# Q_ww (total second derivative of Q(w,Y)=P(w,Y,u1(w,Y),u2(w,Y)) wrt w, at fixed Y -- wait we need wrt w along w only first)
Q_ww = (P_ww + 2*(P_w_u1*u1_w + P_w_u2*u2_w)
        + P_u1u1*u1_w**2 + 2*P_u1u2*u1_w*u2_w + P_u2u2*u2_w**2
        + P_u1*u1_ww + P_u2*u2_ww)

Q_wY = (P_wY + P_w_u1*u1_Y + P_w_u2*u2_Y + P_Y_u1*u1_w + P_Y_u2*u2_w
        + P_u1u1*u1_w*u1_Y + P_u1u2*(u1_w*u2_Y+u2_w*u1_Y) + P_u2u2*u2_w*u2_Y
        + P_u1*u1_wY + P_u2*u2_wY)

Q_YY = (P_YY + 2*(P_Y_u1*u1_Y + P_Y_u2*u2_Y)
        + P_u1u1*u1_Y**2 + 2*P_u1u2*u1_Y*u2_Y + P_u2u2*u2_Y**2
        + P_u1*u1_YY + P_u2*u2_YY)

Q_Y = P_Y + P_u1*u1_Y + P_u2*u2_Y

print("count_ops: Q_ww=",sp.count_ops(Q_ww)," Q_wY=",sp.count_ops(Q_wY)," Q_YY=",sp.count_ops(Q_YY)," Q_Y=",sp.count_ops(Q_Y))

with open('/tmp/Q_second_partials.pkl','wb') as f:
    pickle.dump({'Q_ww':Q_ww,'Q_wY':Q_wY,'Q_YY':Q_YY,'Q_Y':Q_Y}, f)
print("saved")
