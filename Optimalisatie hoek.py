import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pickle
import math

# -----------------------------
# 1. Parameters (SI-eenheden)
# -----------------------------
rho = 7800 #massadichtheid staal
d = 0 #Binnendiameter van het projectiel
D = 16.28 #Buitendiameter van het projectiel
l = 20 #lengte van het projectiel
mu = 1.7 #inschatting
m = rho*(((D/1000)**2-(d/1000)**2)*math.pi/4*l/1000)
L_sp = 90 # Lengte van de spoel
L_ex = 20 #Lengte van de loop na de spoel
end = (L_sp + L_ex - l)/1000 #Positie wanneer de simulatie moet eindigen

A = ((D/1000)**2-(d/1000)**2)*math.pi/4
C_aero = 0.9
h = 1

N=300

dt = 2e-6
n=N/300
#R = 1.804280000757309900
R = 0.68319     # ohm
#L = 2.58947e-3     # H
L = 0.3090207e-3
C = 10e-3             # F
V0 = 62         # V

# -----------------------------
# 2. RLC-simulatie
# -----------------------------
def rlc_discharge(t, x):
    i, di_dt = x
    d2i_dt2 = - (R / L) * di_dt - (1 / (L * C)) * i
    return [di_dt, d2i_dt2]

# Inladen van interpolator
#with open(r"D:\MP\onder gedempt\Data spoel n=300\interpolator_n=300.pkl", "rb") as f:
with open(r"D:\MP\Vol projectiel\Data spoel n=300\interpolator_n=300.pkl", "rb") as f:
#with open(r"D:\MP\Lang\Data spoel n=600\interpolator_n=600.pkl", "rb") as f:
#with open(r"D:\MP\Hol\Data spoel n=600\interpolator_n=600.pkl", "rb") as f:
#with open(r"D:\MP\N var\N=750\interpolator_n=750.pkl", "rb") as f: 
    interpolator = pickle.load(f)

# Bereken de eindpositie voor verschillende hoeken
hoeken = np.linspace(0, 90, 90)  # Hoeken tussen 0° en 45°
beginpositie =-0.00922  # Vaste startpositie in meters
eind_posities = []

for hoek in hoeken:
    print(hoek/90*100,"%")
    hoek_rad = math.radians(hoek)
    F_n = mu * math.cos(hoek_rad) * m * 9.81
    F_g = math.sin(hoek_rad) * m * 9.81
    
    v0 = 0
    i0 = 0
    di_dt0 = V0 / L
    x0_rlc = [i0, di_dt0]

    t_eval = np.arange(0, 0.2, dt)
    sol_rlc = solve_ivp(rlc_discharge, (t_eval[0], t_eval[-1]), x0_rlc, t_eval=t_eval)
    stromen = sol_rlc.y[0]

    posities = np.zeros_like(t_eval)
    snelheden = np.zeros_like(t_eval)
    krachten_coil = np.zeros_like(t_eval)

    posities[0] = beginpositie
    snelheden[0] = v0

    for k in range(1, len(t_eval)):
        if stromen[k] < 1 and stromen[k] > -1:
            F_coil = 0
        else:
            F_coil = interpolator(abs(stromen[k]), posities[k-1] * 1000)

        F_air = 0.5 * 1.29 * A * C_aero * snelheden[k] ** 2
        F_netto = F_coil - F_n - F_g - F_air

        if F_netto < 0 and snelheden[k-1] == 0:
            F_netto = 0

        a = F_netto / m
        snelheden[k] = snelheden[k-1] + a * dt
        posities[k] = posities[k-1] + snelheden[k] * dt
        krachten_coil[k] = F_netto

        if posities[k] >= end:
            posities = posities[:k+1]
            snelheden = snelheden[:k+1]
            krachten_coil = krachten_coil[:k+1]
            stromen = stromen[:k+1]
            t_eval = t_eval[:k+1]
            break
    
    laatste_snelheid = snelheden[-1]
    v_x = laatste_snelheid * math.cos(hoek_rad)
    v_y = laatste_snelheid * math.sin(hoek_rad)
    t_val = (v_y + math.sqrt(v_y**2 + 2 * 9.81 * h)) / 9.81
    eind_pos = v_x * t_val
    eind_posities.append(eind_pos)

# Plot van eindpositie als functie van hoek
plt.figure(figsize=(8, 5))
plt.plot(hoeken, eind_posities, linestyle='-')
plt.xlabel("Hoek (graden)")
plt.ylabel("Eindpositie (m)")
plt.title("Eindpositie als functie van hoek")
plt.grid()
plt.show()

# Bepaal de optimale hoek
beste_hoek = hoeken[np.argmax(eind_posities)]
max_eindpositie = max(eind_posities)
print(f"Beste hoek: {beste_hoek:.2f} graden")
print(f"Maximale eindpositie: {max_eindpositie:.6f} m")
# Eigen data (voorbeeldpunten)
eig_hoeken = [20, 38.43, 41.46, 42.47, 45, 60]                  # graden
eig_afstanden = [1.2, 2, 2.0, 2.1, 1.93, 1.6]            # meters

print("\nVergelijking van eigen hoeken:")
for hoek_eig, afstand_eig in zip(eig_hoeken, eig_afstanden):
    # Zoek de dichtstbijzijnde hoek in de simulatie
    idx = (np.abs(hoeken - hoek_eig)).argmin()
    sim_afstand = eind_posities[idx]
    print(f"Hoek: {hoek_eig:.2f}° | Gemeten: {afstand_eig:.3f} m | Simulatie: {sim_afstand:.3f} m")
    afwijking = (afstand_eig - sim_afstand) / sim_afstand * 100
    print(f"Procentuele afwijking: {afwijking:.2f}%")

# Plot van eindpositie als functie van hoek
plt.figure(figsize=(8, 5))
plt.plot(hoeken, eind_posities, linestyle='-', label='Simulatie')
#plt.plot(eig_hoeken, eig_afstanden, 'x', markersize=5, label='Metingen', color='red')  # Kruisjes
plt.xlabel("Hoek (graden)")
plt.ylabel("Eindpositie (m)")
plt.title("Eindpositie als functie van hoek")
plt.grid()
plt.legend()
plt.show()
