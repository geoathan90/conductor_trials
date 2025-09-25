from sympy import symbols, Eq, solve, re, N
from math import isfinite
import sys

# ---- Input (edit accordingly) ----
a = 290.00          # span (m)
dh = -95.12         # elevation difference h_R - h_L (m). Use 0 for level supports.
w = 1.823           # unit weight (N/m) per HORIZONTAL length
S = 5.47e-4         # area (m^2)
E = 6.18e9          # Young's modulus (Pa)
e = 1.935e-5        # thermal expansion (1/°C)
T1 = 10.0           # initial temp (°C)
T2 = 20.0          # new temp (°C)
H1 = 9.80665 * 2585.0      # initial horizontal tension (N) (from your measured sag)


# ---- Constants ----
ac = e*S*E*(T2-T1)-H1+a**2*S*E*w**2/24/H1**2
bc = a**2*w**2/24
cc = e*(T2-T1)*S*E*w**2*a**2/24 - H1*w**2*a**2/24 - S*E**w**2*a**2/24

# Cubic
H2 = symbols('H2', real=True)
eq = Eq(H2**3 + ac*H2**2 + bc*H2 + cc, 0)

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
# Sag relative to each support (parabolic)
fL = (w / (2.0 * H2_sol)) * x0 * (S - x0)
fR = fL  # symmetric formula; clearances differ due to support heights

# print(f"Sag at mid (level supports)   f = {f_mid:.3f} m")
# print(f"Lowest point from left        x0 = {x0:.3f} m  (useful when dh ≠ 0)")
# print(f"Sag rel. to left support     fL = {fL:.3f} m")
# print(f"Sag rel. to right support    fR = {fR:.3f} m")
