import numpy as np
from scipy.integrate import solve_ivp
import pickle
import math
import matplotlib.pyplot as plt

# -----------------------------
# Parameters (SI units)
# -----------------------------
rho = 7800                   # kg/m^3
D = 16                       # mm
d = 14                       # mm
l = 20                       # mm
mu = 4                       # wrijving

# massa van de spoel
area = ((D / 1000) ** 2 - (d / 1000) ** 2) * math.pi / 4
m = rho * area * (l / 1000)

L_sp = 90                    # mm
L_ex = 30                    # mm
end = (L_sp + L_ex) / 1000   # m

hoek_deg = 38.43             # graden
pos_0 = -0.013342            # m
hoek = math.radians(hoek_deg)

A = area                     # m^2
C_aero = 0.9
h = 1                        # m
g = 9.81

F_n = mu * math.cos(hoek) * m * g
F_g = math.sin(hoek) * m * g

# RLC circuit
R = 0.68319
L = 0.3090207e-3
C = 10e-3
V0 = 62

dt = 2e-6
t_total = 0.0175
t_eval = np.arange(0, t_total, dt)

# Interpolator
def load_interpolator(path):
    with open(path, "rb") as f:
        return pickle.load(f)

interpolator = load_interpolator(r"D:\MP\onder gedempt\Data spoel n=300\interpolator_n=300.pkl")

# -----------------------------
# Simulatie
# -----------------------------
def simulate_full(t_on):
    global stromen, krachten_coil, posities, snelheden  # nodig voor plotten

    # RLC ontlading
    def rlc_discharge(t, x):
        i, di_dt = x
        return [di_dt, - (R / L) * di_dt - (1 / (L * C)) * i]

    x0_rlc = [0, V0 / L]
    sol = solve_ivp(rlc_discharge, (0, t_total), x0_rlc, t_eval=t_eval)
    stromen = sol.y[0].copy()
    stromen[t_eval > t_on] = 0

    # Initialisatie
    posities = np.zeros_like(t_eval)
    snelheden = np.zeros_like(t_eval)
    krachten_coil = np.zeros_like(t_eval)
    krachten_netto = np.zeros_like(t_eval)

    posities[0] = pos_0

    for k in range(1, len(t_eval)):
        stroom = abs(stromen[k])
        positie = posities[k - 1]

        if stroom < 1:
            F_coil = 0
        else:
            F_coil = interpolator(stroom, positie * 1000)

        krachten_coil[k] = F_coil

        F_air = 0.5 * 1.29 * A * C_aero * snelheden[k - 1] ** 2
        direction = -1 if snelheden[k - 1] > 0 else 1

        F_tegen = F_n + F_g + F_air
        F_netto = F_coil + direction * F_tegen

        if F_netto < 0 and snelheden[k - 1] == 0:
            F_netto = 0

        krachten_netto[k] = F_netto

        snelheden[k] = snelheden[k - 1] + (F_netto / m) * dt
        posities[k] = positie + snelheden[k] * dt

        if posities[k] >= end:
            v_exit = snelheden[k]
            break
    else:
        v_exit = snelheden[-1]

    v_x = v_exit * math.cos(hoek)
    v_y = v_exit * math.sin(hoek)
    t_vlucht = (v_y + math.sqrt(v_y ** 2 + 2 * g * h)) / g
    afstand = v_x * t_vlucht

    return afstand, v_exit, np.max(np.abs(stromen)), np.max(krachten_netto)

# Run en plot
if __name__ == "__main__":
    afstand, v_exit, i_max, F_max = simulate_full(t_on=0.4)
    print(f"Exit velocity: {v_exit:.6f} m/s")
    print(f"Horizontal distance: {afstand:.6f} m")
    print(f"Max current: {i_max:.2f} A")
    print(f"Max net force: {F_max:.2f} N")

    # Figuren
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(4, 1, figsize=(10, 8), sharex=True)

    axs[0].plot(t_eval, stromen, 'b', label="Stroom (A)")
    axs[0].set_ylabel("Stroom (A)")
    axs[0].set_title("Stroom als functie van de tijd")
    axs[0].grid(True)
    axs[0].legend()

    axs[1].plot(t_eval, krachten_coil, 'r', label="Kracht (N)")
    axs[1].set_ylabel("Kracht (N)")
    axs[1].set_title("Kracht als functie van de tijd")
    axs[1].grid(True)
    axs[1].legend()

    axs[2].plot(t_eval, posities, 'g', label="Positie (m)")
    axs[2].set_ylabel("Positie (m)")
    axs[2].set_title("Positie als functie van de tijd")
    axs[2].grid(True)
    axs[2].legend()

    axs[3].plot(t_eval, snelheden, color="purple", label="Snelheid (m/s)")
    axs[3].set_xlabel("Tijd (s)")
    axs[3].set_ylabel("Snelheid (m/s)")
    axs[3].set_title("Snelheid als functie van de tijd")
    axs[3].grid(True)
    axs[3].legend()

    plt.tight_layout()
    plt.show()

    # Kracht vs positie
    plt.figure(figsize=(8, 5))
    plt.plot(posities, krachten_coil, color="orange", label="Kracht (N) vs Positie (m)")
    plt.xlabel("Positie (m)")
    plt.ylabel("Kracht (N)")
    plt.title("Kracht als functie van de positie")
    plt.grid(True)
    plt.legend()
    plt.show()
