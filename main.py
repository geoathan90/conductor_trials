# sympy_cubic_stringing.py
# Minimal SymPy solve for H2, then sag. Parabolic model.

from sympy import symbols, Eq, solve, re, N
from math import isfinite, sinh, asinh
import sys

# ---- Input (edit these) ----
S = 716.72         # span (m)
dh = 37.54        # elevation difference h_R - h_L (m). Use 0 for level supports.
w = 1.823          # unit weight (N/m) per HORIZONTAL length
A = 5.47e-4        # area (m^2)
E = 6.18e9         # Young's modulus (Pa)
alpha = 1.935e-5     # thermal expansion (1/°C)
T1 = 50.0          # initial temp (°C)
T2 = -18.0         # new temp (°C)
H1 = 9.80665 * 2585.0
H1_old = H1      # initial horizontal tension (N) (from your measured sag)
#H1 = float(sys.argv[1])*9.80665

# ---- Constants ----
b = 1.0 / (E * A)
k = (w**2 * S**2) / 24.0
dT = T2 - T1

# Cubic: b*H2^3 + (alpha*dT - b*H1 + k/H1^2)*H2^2 - k = 0
H2 = symbols('H2', real=True)
c2 = alpha * dT - b * H1 + (k / (H1**2))
eq = Eq(b*H2**3 + c2*H2**2 - k, 0)

roots = solve(eq, H2)  # symbolic roots
# Pick the physically valid root: real, positive, finite
candidates = []
for r in roots:
    rv = complex(N(r))
    if abs(rv.imag) < 1e-9 and rv.real > 0 and isfinite(rv.real):
        candidates.append(rv.real)

if not candidates:
    raise RuntimeError("No positive real root for H2. Check inputs.")

H2_sol = max(candidates)  # tensions are typically the largest positive real root
print(f"H1 = {H1*0.102:.3f} kp")
print(f"H2 = {H2_sol*0.102:.3f} kp")

###############################################


# ---- Sag (parabolic) ----
# If supports are level: f = w S^2 / (8 H2)
f_mid = w * S**2 / (8.0 * H2_sol)

# If supports are at different heights, the lowest point shifts:
# x0 = S/2 + (H2*dh)/(w*S)
x0 = S/2.0 + (H2_sol * dh) / (w * S)
x0_old = S/2.0 + (H1_old * dh) / (w * S)

x0 = -H2_sol/w*asinh(dh/2/(H2_sol/w)/sinh(S/2/(H2_sol/w)))+S/2
x0_old = -H1_old/w*asinh(dh/2/(H1_old/w)/sinh(S/2/(H1_old/w)))+S/2

# Sag relative to each support (parabolic)
fL = (w / (2.0 * H2_sol)) * x0 * (S - x0)
fR = fL  # symmetric formula; clearances differ due to support heights

# print(f"Sag at mid (level supports)   f = {f_mid:.3f} m")
# print(f"Lowest point from left for H1       x0_old = {x0:.3f} m  ")
# print(f"Lowest point from left for H2       x0 = {x0:.3f} m  ") #(useful when dh ≠ 0)
# print(f"Sag rel. to left support     fL = {fL:.3f} m")
# print(f"Sag rel. to right support    fR = {fR:.3f} m")
