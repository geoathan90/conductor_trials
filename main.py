
from sympy import symbols, Eq, solve, re, N
from math import isfinite, sinh, asinh
import sys

# ---- Input ----
S =  50           # span (m)      290 316.72
dh =  0             # elevation difference h_R - h_L (m)     -95.12  -37.54
w = 1.823           #
w1 = w              # kg/m
w2 = w              #    
A = 5.47e-4         # area (m^2)
E = 5.132e9          # Young's modulus (Pa)  6.184e9
alpha = 1.935e-5    # thermal expansion (1m/°C)
T1 = 0          # initial temp (°C)
T2 = 40          # new temp (°C)
H1 = 442           # initial horizontal tension (kg)   ex:2585 για τις μαλακίες τους
H1_old = H1         #  archive
#H1 = float(sys.argv[1])

S1 = 377.52 #pyrgos 22
dh1 = -85.42
S2 = 369.70 #pyrgos 25
dh2 = -70.70

### note: νομίζω έχουν υπολογίσει τα βάρη στους -10 με γυμνό αγωγό 

def solve_for_H2(S,H1,T1,T2,w1=w,w2=w):

    # Coefficients
    c1 = alpha * A * E * (T2 - T1) - H1 + (w1**2 * A * E * S**2) / (24.0 * H1**2)
    c2 = S**2 * w2**2 / 24.0
    c3 = alpha * A * E * (T2 - T1) * (S**2 * w2**2) / 24.0 - H1 * (S**2 * w2**2) / 24.0 - (w2**2 * A * E * S**2) / 24.0

    # Cubic: x^3 + c1*x^2 + c2*x + c3 = 0
    H2 = symbols('H2')  
    eq = Eq(H2**3 + c1*H2**2 + c2*H2 + c3, 0)   

    roots = solve(eq, H2)  
    #print(f"Found roots: {roots}")
    
    # pick valid root: real, positive, finite
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

def Th_from_sag(sag, S, w):
    return w * S**2 / (8.0 * sag)    

def distance_lowest_point_r(S,dh,H=H1,w=w):
    return -H/w*asinh(dh/2/(H/w)/sinh(S/2/(H/w)))+S/2

def distance_lowest_point_l(S,dh,H=H1,w=w):
    dh = -dh
    return -H/w*asinh(dh/2/(H/w)/sinh(S/2/(H/w)))+S/2

def monopleyro_right(S,dh,H=H1,w=w):
    monopleyro = -H/w*asinh(dh/2/(H/w)/sinh(S/2/(H/w)))-S/2
    if monopleyro <= 0: print("δεν υπάρχει μονόπλευρο φορτίο") 
    else: return monopleyro

def monopleyro_left(S,dh,H=H1,w=w):
    dh = -dh
    monopleyro = -H/w*asinh(dh/2/(H/w)/sinh(S/2/(H/w)))-S/2
    if monopleyro <= 0: print("δεν υπάρχει μονόπλευρο φορτίο") 
    else: return monopleyro

def synoliko_katakoryfo(S_l,dh_l,H_l,S_r,dh_r,H_r):
    return distance_lowest_point_r(S_r,dh_r,H_r) + distance_lowest_point_l(S_l,dh_l,H_l)

#print(f"Sag at mid (level supports)   f = {f_mid:.3f} m")
#print(f"x0_old = {x0_old:.3f} m  ")
#print(f"x0 = {x0:.3f} m  ") #(useful when dh ≠ 0)

# examples (290,-95.12,2585), (316.72,37.54,2585), (377.52,-85.42,2585)