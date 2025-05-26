# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 19:30:02 2025
@author: Tibo
"""

import numpy as np
from scipy.integrate import solve_ivp
import pickle
import math

# -----------------------------
# 1. Parameters (SI-eenheden)
# -----------------------------
rho = 7800                   # dichtheid (kg/m^3)
d = 14                       # mm, binnendiameter spoel
D = 16                       # mm, buitendiameter spoel
l = 20                       # mm, spoellengte
mu = 4                       # wrijvingscoëfficiënt

# massa berekenen: volume maal dichtheid
m = rho * (((D/1000)**2 - (d/1000)**2) * math.pi/4 * l/1000)

L_sp = 90                    # mm
L_ex = 30                    # mm
end = (L_sp + L_ex) / 1000   # Eindpositie in meters

hoek_deg = 41.46             # hoek in graden, als launchhoek
hoek = math.radians(hoek_deg)  # omrekenen naar radialen

A = (((D/1000)**2 - (d/1000)**2) * math.pi/4)  # dwarsdoorsnede spoel
C_aero = 0.9                # aerodynamische factor

h = 1                       # beginhoogte in meter

# Natuurkundige krachten (zoals wrijvingskracht op de spoel)
F_n = mu * math.cos(hoek) * m * 9.81
F_g = math.sin(hoek) * m * 9.81

# Parameters voor het RLC-circuit
R = 0.68319         # ohm
L = 0.3090207e-3
C = 10e-3                   # F
V0 = 62                     # V
dt = 2e-6                   # tijdstap in s

# Laad de interpolator voor de spoelkrachten
with open(r"D:\MP\onder gedempt\Data spoel n=300\interpolator_n=300.pkl", "rb") as f:
#with open(r"D:\MP\Hol\Data spoel n=600\interpolator_n=600.pkl", "rb") as f:
#with open(r"D:\MP\Hol\Data spoel n=600\interpolator_n=600.pkl", "rb") as f:
    interpolator = pickle.load(f)


def simulate_trajectory(t_on, hoek=hoek, t_total=0.2):
    """Simuleert de spoelaandrijving + projectieltraject voor een gegeven t_on."""
    t_eval = np.arange(0, t_total, dt)

    # RLC-simulatie
    def rlc_discharge(t, x):
        i, di_dt = x
        return [di_dt, - (R / L) * di_dt - (1 / (L * C)) * i]

    x0_rlc = [0, V0 / L]
    sol_rlc = solve_ivp(rlc_discharge, (0, t_total), x0_rlc, t_eval=t_eval)
    stromen = sol_rlc.y[0].copy()
    stromen[t_eval > t_on] = 0

    # Mechanische simulatie
    pos = np.zeros_like(t_eval)
    vel = np.zeros_like(t_eval)
    krachten = np.zeros_like(t_eval)

    pos[0] = -0.02
    for k in range(1, len(t_eval)):
        # coil-kracht
        if abs(stromen[k]) < 1:
            F_coil = 0
        else:
            F_coil = interpolator(abs(stromen[k]), pos[k-1]*1000)

        # luchtweerstand
        F_air = 0.5 * 1.29 * A * C_aero * vel[k-1]**2
        fact = -1 if vel[k-1] > 0 else 1

        F_tegen = F_n + F_g + F_air
        F_netto = F_coil + fact * F_tegen
        if F_netto < 0 and vel[k-1] == 0:
            F_netto = 0

        a = F_netto / m
        vel[k] = vel[k-1] + a * dt
        pos[k] = pos[k-1] + vel[k] * dt
        krachten[k] = F_netto

        if pos[k] >= end:
            pos = pos[:k+1]
            vel = vel[:k+1]
            t_eval = t_eval[:k+1]
            break

    v_final = vel[-1]
    v_x = v_final * math.cos(hoek)
    v_y = v_final * math.sin(hoek)
    t_val = (v_y + math.sqrt(v_y**2 + 2 * 9.81 * h)) / 9.81
    delta_x = v_x * t_val

    return delta_x


def bracket_for_target(target, t_min=0, t_max=0.2, n=30):
    """
    Zoekt binnen [t_min, t_max] naar een subinterval [t_i, t_{i+1}]
    waar simulate_trajectory(t_i) <= target <= simulate_trajectory(t_{i+1}).
    """
    ts = np.linspace(t_min, t_max, n)
    xs = [simulate_trajectory(t) for t in ts]
    for i in range(1, len(ts)):
        if xs[i-1] <= target <= xs[i]:
            print(f"Bracket: t_on ∈ [{ts[i-1]:.6f}, {ts[i]:.6f}] "
                  f"met afstanden [{xs[i-1]:.4f}, {xs[i]:.4f}] m")
            return ts[i-1], ts[i]
    raise ValueError("Doelafstand buiten bereik.")


def find_t_on_for_target(target_distance, tol=1e-5):
    """
    Berekent t_on m.b.v. bisectie, na eerst automatisch bracket te vinden.
    """
    t_on_min, t_on_max = bracket_for_target(target_distance)
    

    while t_on_max - t_on_min > tol:
        t_on_mid = 0.5 * (t_on_min + t_on_max)
        pos_mid = simulate_trajectory(t_on_mid)
        print(f" t_on={t_on_mid:.6f} → Δx={pos_mid:.6f} m")
        if pos_mid < target_distance:
            t_on_min = t_on_mid
        else:
            t_on_max = t_on_mid

    return 0.5 * (t_on_min + t_on_max)


if __name__ == "__main__":
    try:
        target_distance = float(input("Voer gewenste horizontale afstand in (m): "))
    except ValueError:
        print("Ongeldige invoer, voer een getal in meters in.")
        exit(1)

    print(f"Zoeken naar t_on voor target = {target_distance:.2f} m …")
    t_on_target = find_t_on_for_target(target_distance)
    eind_afstand = simulate_trajectory(t_on_target)

    print(f"\nBenodigde aan-tijd (t_on): {t_on_target*10e5:.0f} µs")
    print(f"Behaalde afstand: {eind_afstand:.2f} m")


