
from sympy import symbols, Eq, solve, re, N
from math import isfinite, sinh, asinh
import sys

# ---- Input ----
S =  377.52           # span (m)      290 316.72
dh =  -85.42             # elevation difference h_R - h_L (m)     -95.12  -37.54
w = 1.823           # kg/m
A = 5.47e-4         # area (m^2)
E = 6.18e9          # Young's modulus (Pa)
alpha = 1.935e-5    # thermal expansion (1m/°C)
T1 = 50          # initial temp (°C)
T2 = -17.7778          # new temp (°C)
H1 = 2585           # initial horizontal tension (kg)   ex:2585 για τις μαλακίες τους
H1_old = H1         #  archive
#H1 = float(sys.argv[1])

def solve_for_H2(S,H1,T2,T1):
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
        raise RuntimeError("Δεν βρέθηκε θετική ρίζα, τσέκαρε τα input.")

    H2_sol = max(candidates)  # tensions are typically the largest positive real root

    return H2_sol


def sag(S, H):
    return w * S**2 / (8.0 * H)


def distance_lowest_point(S,dh,H=H1):
    return -H/w*asinh(dh/2/(H/w)/sinh(S/2/(H/w)))+S/2

def monopleyro_right(S,dh,H=H1):
    monopleyro = -H/w*asinh(dh/2/(H/w)/sinh(S/2/(H/w)))-S/2
    if monopleyro <= 0: print("δεν υπάρχει μονόπλευρο φορτίο") 
    else: return monopleyro

def monopleyro_left(S,dh,H=H1):
    dh = -dh
    monopleyro = -H/w*asinh(dh/2/(H/w)/sinh(S/2/(H/w)))-S/2
    if monopleyro <= 0: print("δεν υπάρχει μονόπλευρο φορτίο") 
    else: return monopleyro

def synoliko_katakoryfo(S_r,dh_r,H_r,S_l,dh_l,H_l):
    dh_l = -dh_l
    return distance_lowest_point(S_r,dh_r,H_r) + distance_lowest_point(S_l,dh_l,H_l) 

# Sag relative to each support (parabolic)
#fL = (w / (2.0 * H2_sol)) * x0 * (S - x0)
#fR = fL  # symmetric formula; clearances differ due to support heights

#print(f"Sag at mid (level supports)   f = {f_mid:.3f} m")
#print(f"x0_old = {x0_old:.3f} m  ")
#print(f"x0 = {x0:.3f} m  ") #(useful when dh ≠ 0)
# print(f"Sag rel. to left support     fL = {fL:.3f} m")
# print(f"Sag rel. to right support    fR = {fR:.3f} m")

# examples (290,-95.12,2585), (316.72,37.54,2585), (377.52,-85.42,2585)