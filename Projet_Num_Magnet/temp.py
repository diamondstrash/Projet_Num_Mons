import matplotlib.pyplot as plt
import numpy as np
import cmath
import math

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

# Calculer l'argument de chaque élément de la liste 'Impedance' en radians
argument_radians = [cmath.phase(z) for z in Impedance]

# Convertir l'argument en degrés
phi = [math.degrees(arg) for arg in argument_radians]

# Create the figure and the primary y-axis
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot resistance (Ohms) on the primary y-axis
ax1.plot(frequencies, R_ohm, label="Resistance (Ω)", color="blue")
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("Resistance (Ω)", color="blue")
ax1.tick_params(axis='y', labelcolor="blue")

# Create a secondary y-axis for power (W)
ax2 = ax1.twinx()
ax2.plot(frequencies, Pcoil, label="Power (W)", color="red")
ax2.set_ylabel("Power (W)", color="red")
ax2.tick_params(axis='y', labelcolor="red")

# Create another secondary y-axis for energy (J)
ax3 = ax1.twinx()
ax3.spines['right'].set_position(('outward', 60))  # Move the axis to the right
ax3.plot(frequencies, Wm, label="Energy (J)", color="green")
ax3.set_ylabel("Energy (J)", color="green")
ax3.tick_params(axis='y', labelcolor="green")

# Create a fourth y-axis for inductance (H)
ax4 = ax1.twinx()
ax4.spines['right'].set_position(('outward', 120))  # Move the axis to the right
ax4.plot(frequencies, L, label="Inductance (H)", color="purple")
ax4.set_ylabel("Inductance (H)", color="purple")
ax4.tick_params(axis='y', labelcolor="purple")

# Create a fifth y-axis for core power (Pcore)
ax5 = ax1.twinx()
ax5.spines['right'].set_position(('outward', 180))  # Move the axis further to the right
ax5.plot(frequencies, Pcore, label="Core Power (W)", color="orange")
ax5.set_ylabel("Core Power (W)", color="orange")
ax5.tick_params(axis='y', labelcolor="orange")

# Add a legend
fig.legend(loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=5)

# Adjust layout for better clarity
fig.tight_layout()

# Set the title
plt.title("Multi-Unit Graph with Multiple Y-Axes")

# Show the plot
plt.show()