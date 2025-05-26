import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pickle
import math

# -----------------------------
# 1. Parameters (SI-eenheden)
# -----------------------------
rho = 7800
D = 16      # mm, buitendiameter
d = 14      # mm, binnendiameter
l = 20      # mm, spoellengte
mu = 4    # wrijvingscoëfficiënt

# massa van de spoel
area = ((D/1000)**2 - (d/1000)**2) * math.pi/4
m = rho * area * (l/1000)

L_sp = 90   # mm, spoelsteek
L_ex = 20   # mm, extra slag
end = (L_sp + L_ex - l)/1000  # m, eindpositie

hoek = 41.46            # launch hoek in graden
hoek = math.radians(hoek)

A = area            # m^2, dwarsdoorsnede spoel
C_aero = 0.9        # aerodynamische factor
h = 1               # m, beginhoogte

g = 9.81           # m/s^2, zwaartekracht
F_n = mu * math.cos(hoek) * m * g
F_g = math.sin(hoek) * m * g

# RLC-circuit parameters
R = 0.68319    # ohm
L = 0.3090207e-3     # H
C_circ = 10e-3     # F (capaciteit)
V0 = 62            # V

dt = 2e-6          # s, tijdstap

# laad interpolator
def load_interpolator(path):
    with open(path, "rb") as f:
        return pickle.load(f)

interpolator = load_interpolator(r"D:\MP\onder gedempt\Data spoel n=300\interpolator_n=300.pkl")

# -----------------------------
# 2. RLC-simulatie
# -----------------------------
def rlc_discharge(t, x):
    i, di_dt = x
    return [di_dt, - (R/L) * di_dt - (1/(L*C_circ)) * i]

# -----------------------------
# 3. Startpositie-scan
# -----------------------------
start_positions = np.linspace(-0.02, 0.0075, 100)
eind_posities = []

t_eval_full = np.arange(0, 0.2, dt)

for x0 in start_positions:
    # RLC
    sol = solve_ivp(rlc_discharge,
                    (t_eval_full[0], t_eval_full[-1]),
                    [0, V0/L],
                    t_eval=t_eval_full)
    stromen = sol.y[0]

    # mechanica
    pos = np.zeros_like(t_eval_full)
    vel = np.zeros_like(t_eval_full)
    pos[0] = x0

    for k in range(1, len(t_eval_full)):
        # spoelkracht
        if abs(stromen[k]) < 1:
            F_coil = 0
        else:
            F_coil = interpolator(abs(stromen[k]), pos[k-1]*1000)
        # break bij NaN
        if np.isnan(F_coil):
            pos = pos[:k]
            vel = vel[:k]
            stromen = stromen[:k]
            t_eval = t_eval_full[:k]
            break
        
        # luchtweerstand
        F_air = 0.5 * 1.29 * A * C_aero * vel[k-1]**2
        F_netto = F_coil - (F_n + F_g + F_air)
        if F_netto < 0 and vel[k-1] == 0:
            F_netto = 0

        # update
        vel[k] = vel[k-1] + F_netto/m * dt
        pos[k] = pos[k-1] + vel[k] * dt

        # stop bij eindpositie
        if pos[k] >= end:
            pos = pos[:k+1]
            vel = vel[:k+1]
            stromen = stromen[:k+1]
            t_eval = t_eval_full[:k+1]
            break

    # eindsnelheid
    v_exit = vel[-1]
    # projectiel
    v_x = v_exit * math.cos(hoek)
    v_y = v_exit * math.sin(hoek)
    t_val = (v_y + math.sqrt(v_y**2 + 2*g*h)) / g
    eind_posities.append(v_x * t_val)

# -----------------------------
# 4. Plot en resultaat
# -----------------------------
plt.figure(figsize=(8,5))
plt.plot(start_positions, eind_posities)
plt.xlabel("Startpositie (m)")
plt.ylabel("Eindpositie (m)")
plt.title("Eindpositie vs startpositie")
plt.grid(True)
plt.show()

beste = start_positions[np.argmax(eind_posities)]
print(f"Beste startpositie: {beste:.6f} m")
print(f"Maximale eindpositie: {max(eind_posities):.6f} m")
# Voeg eigen data toe (voorbeeldpunten)


eigen_startposities = [-0.02, -0.01171, -0.005, 0, 0.005]   # in meters
eigen_eindposities = [1.42, 4.2, 2.77, 1.46, 0.55]            # in meters

# Plot simulatiecurve
plt.plot(start_positions, eind_posities, label='Simulatie', linestyle='-', color='blue')

# Plot eigen data met kruisjes
plt.plot(eigen_startposities, eigen_eindposities, 'x', markersize=5, label='Meetpunten', color='red')

# Labels en legende
plt.xlabel("Startpositie (m)")
plt.ylabel("Eindpositie (m)")
plt.title("Eindpositie vs startpositie")
plt.grid(True)
plt.legend()

from scipy.interpolate import interp1d

# Interpoleer de simulatiecurve
simulatie_interp = interp1d(start_positions, eind_posities, kind='linear', fill_value="extrapolate")

print("\nVergelijking van eigen startposities:")
for x_meet, y_meet in zip(eigen_startposities, eigen_eindposities):
    y_sim = simulatie_interp(x_meet)
    afwijking = 100 * (y_meet - y_sim) / y_sim
    print(f"Startpositie: {x_meet:.4f} m | Simulatie: {y_sim:.3f} m | Meting: {y_meet:.3f} m → Afwijking: {afwijking:.2f} %")
