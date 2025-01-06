import matplotlib.pyplot as plt
import cmath
import math
import numpy as np
import os
from scipy.optimize import curve_fit

# Données
S_coil = 1e-5
N_turns = 125
sigma_Cu = 5.96e7
rb = 19e-3
b = 1e-3
l = 10e-3

R_coil = ((N_turns ** 2)*(2*math.pi)*(rb+b/2))/(sigma_Cu*S_coil)

frequencies = [25, 50, 100, 500, 1000]
omega = [2 * math.pi * f for f in frequencies]

# Puissance et énergies
Pcore = np.array([0.062355, 0.047348, 0.0351, 0.016393, 0.011538]) 
Pcoil = np.array([0.016848, 0.008969, 0.004674, 9.685034e-4, 4.811176e-4])
Wm = np.array([4.134096e-4, 1.563597e-4, 5.799992e-5, 5.52802e-6, 1.982118e-6])

# Liste de nombres complexes
densite_courant = [9.900357e5-8.117285e5j, 7.039679e5-6.140233e5j, 4.971708e5 - 4.55531e5j, 2.170225e5 - 2.170852e5j, 1.502451e5 - 1.556755e5j]

densite_courant_norm = [abs(z) for z in densite_courant]

courant = [(S_coil / N_turns ) * z for z in densite_courant]

# Diviser 2 par chaque élément de la liste
Impedance = [2 / z for z in courant]

# Parties réelles et imaginaires
R_ohm = [z.real for z in Impedance]
X_ohm = [z.imag for z in Impedance]

# Calcul de L
L = [x / o for x, o in zip(X_ohm, omega)] 

R_core_bis = np.zeros(len(R_ohm))

for i in range(len(R_ohm)):
    R_core_bis[i] = (R_ohm[i] - R_coil) / (N_turns ** 2)

print(R_core_bis)

# Calculer l'argument de chaque élément de la liste 'Impedance' en radians
argument_radians = [cmath.phase(z) for z in Impedance]

# Convertir l'argument en degrés
phi = [math.degrees(arg) for arg in argument_radians]

# Plotting grouped by unit

# Create a folder to save figures if it doesn't exist
output_folder = "Figure"
os.makedirs(output_folder, exist_ok=True)

# 1. Power-related quantities (Pcore, Pcoil)
plt.figure(figsize=(10, 6))
plt.plot(frequencies, Pcore, label="Pcore (W)", marker="o")
plt.plot(frequencies, Pcoil, label="Pcoil (W)", marker="o")
plt.title("Power as a function of frequency")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power (W)")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_folder, "power_vs_frequency.png"))
plt.close()  # Close the figure to free memory

# 2. Impedance-related quantities (R_ohm, X_ohm, R_core_bis)
plt.figure(figsize=(10, 6))
plt.plot(frequencies, R_ohm, label="R (Ohms)", marker="o")
plt.plot(frequencies, X_ohm, label="X (Ohms)", marker="o")
plt.title("Impedance-related quantities as a function of frequency")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Impedance (Ohms)")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_folder, "impedance_vs_frequency.png"))
plt.close()

# 3. Inductance-related quantities (L)
plt.figure(figsize=(10, 6))
plt.plot(frequencies, L, label="L (H)", marker="o")
plt.title("Inductance as a function of frequency")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Inductance (H)")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_folder, "inductance_vs_frequency.png"))
plt.close()

# 4. Phase-related quantity (phi)
plt.figure(figsize=(10, 6))
plt.plot(frequencies, phi, label=r"$\phi$ (degrees)", marker="o")
plt.title("Phase angle as a function of frequency")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Phase angle (°)")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_folder, "phase_vs_frequency.png"))
plt.close()

# 5. Magnetic energy (Wm)
plt.figure(figsize=(10, 6))
plt.plot(frequencies, Wm, label="Wm (J)", marker="o", color="purple")
plt.title("Magnetic energy as a function of frequency")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnetic energy (J)")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_folder, "magnetic_energy_vs_frequency.png"))
plt.close()

# 6. Additional calculations (L_wm, R_core)
courant_norm = [abs(z) for z in courant]

L_wm = [2 * x / (o ** 2) for x, o in zip(Wm, courant_norm)] #calcul de l'inductance avec Wm et I^2 (norme)

R_core = [2 * x / ((N_turns*o)**2) for x, o in zip(Pcore, courant_norm)] #calcul de R_core avec P_core et I^2