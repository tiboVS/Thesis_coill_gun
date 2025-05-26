import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, cumtrapz

N=300

# Parameters van de kring
R = 0.68319        # ohm
#L = 0.001623275389530026 
L = 0.413965e-3
C = 10e-3  # Capaciteit in Farad

# Systeem van ODE's zonder bron
def rlc_discharge(t, x):
    i, di_dt = x
    d2i_dt2 = - (R / L) * di_dt - (1 / (L * C)) * i
    return [di_dt, d2i_dt2]

# Beginwaarden
V0 = 60  # Beginspanning over condensator in Volt
i0 = 0   # Beginstroom
di_dt0 = V0 / L  # di/dt(0) = V0 / L

x0 = [i0, di_dt0]

# Simulatieperiode
t_span = (0, 0.02)
t_eval = np.linspace(t_span[0], t_span[1], 1000)

# Oplossen
sol = solve_ivp(rlc_discharge, t_span, x0, t_eval=t_eval)

# Bereken spanning over condensator door integratie van stroom
Vc = V0 - (1/C) * cumtrapz(sol.y[0], sol.t, initial=0)

# Gewenste tijdstippen om de stroom te bepalen
t_points = [0.005, 0.01, 0.015]

# Interpoleer de stroom op deze tijdstippen
i_values = np.interp(t_points, sol.t, sol.y[0])

# Print stroom op opgegeven tijdstippen
for t, i in zip(t_points, i_values):
    print(f"Stroom op t = {t:.3f} s: {i:.4f} A")

# Maximumwaarden
i_max = np.max(sol.y[0])
t_i_max = sol.t[np.argmax(sol.y[0])]

Vc_max = np.max(Vc)
t_Vc_max = sol.t[np.argmax(Vc)]

print(f"\nMaximale stroom: {i_max:.4f} A op t = {t_i_max:.6f} s")
print(f"Maximale spanning over condensator: {Vc_max:.4f} V op t = {t_Vc_max:.6f} s")

# Plotten met tweede y-as
fig, ax1 = plt.subplots()

# Eerste y-as (stroom)
color = 'tab:blue'
ax1.set_xlabel("Tijd (s)")
ax1.set_ylabel("Stroom i(t) [A]", color=color)
ax1.plot(sol.t, sol.y[0], color=color, label='Stroom i(t)')
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid(True)

# Tweede y-as (spanning)
ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel("Spanning over condensator Vc(t) [V]", color=color)
ax2.plot(sol.t, Vc, color=color, linestyle='--', label='Vc(t)')
ax2.tick_params(axis='y', labelcolor=color)

plt.title("Ontlading van Serie RLC-kring (Stroom en Condensatorspanning)")
plt.tight_layout()
plt.show()
