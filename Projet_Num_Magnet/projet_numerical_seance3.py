import matplotlib.pyplot as plt
import cmath
import math
import numpy as np

# Données
S_coil = 1e-5
N_turns = 125
sigma_Cu = 5.96e7
rb = 19e-3
b = 1e-3
l = 10e-3

R_coil = ((N_turns ** 2)*(2*math.pi)*(rb+b/2))/(sigma_Cu*S_coil)
print(R_coil)

frequencies = [25, 50, 100, 500, 1000]
omega = [2 * math.pi * f for f in frequencies]

#Puissance et énergies
Pcore = [0.062355, 0.047348, 0.0351, 0.016393, 0.011538]
Pcoil = [0.016848, 0.008969, 0.004674, 9.685034e-4, 4.811176e-4]
Wm = [4.134096e-4, 1.563597e-4, 5.799992e-5, 5.52802e-6, 1.982118e-6]

# Liste de nombres complexes
densite_courant = [9.900357e5-8.117285e5j,7.039679e5-6.140233e5j,4.971708e5 - 4.55531e5j,2.170225e5 - 2.170852e5j, 1.502451e5 - 1.556755e5j]

courant = [(S_coil / N_turns ) * z for z in densite_courant]
courant_core = [-9.900354+8.117285j, -7.039674+6.140232j, -4.971716+4.555305j, -2.170222+2.170844j, -1.502431+1.556735j]

# Diviser 2 par chaque élément de la liste
Impedance = [2 / z for z in courant]

# Parties réelles et imaginaires
R_ohm = [z.real for z in Impedance]
X_ohm = [z.imag for z in Impedance]

# Calcul de L
L = [x / o for x, o in zip(X_ohm, omega)] 

R_core_bis = np.zeros(len(R_ohm))

for i in range(len(R_ohm)):
    R_core_bis[i] = (R_ohm[i] - R_coil)/(N_turns**2)

# Calculer l'argument de chaque élément de la liste 'Impedance' en radians
argument_radians = [cmath.phase(z) for z in Impedance]

# Convertir l'argument en degrés
phi = [math.degrees(arg) for arg in argument_radians]

# Création du graphique
plt.figure(figsize=(10, 6))

# Tracer les courbes
plt.plot(frequencies, Pcore, label="Pcore (W)", marker="o")
plt.plot(frequencies, Wm, label="Wm (J)", marker="o")
plt.plot(frequencies, R_ohm, label="R (ohm)", marker="o")
plt.plot(frequencies, X_ohm, label="X (ohm)", marker="o")
plt.plot(frequencies, R_core_bis, label=r'$R_{core}$ (ohm)', marker="o")
plt.plot(frequencies, L, label="L (H)", marker="o")
plt.plot(frequencies, phi, label="\u03C6 (°)", marker="o")

# Personnalisation
plt.title("Evolution of the different parameters as a function of frequency")
plt.xlabel("Frequency (f)")
plt.ylabel("Values")
plt.legend()
plt.grid(True)

# Afficher le graphique
plt.show()

# Calcul de la norme pour chaque nombre complexe (courant)
courant_norm = [abs(z) for z in courant]
courant_core = [abs(z) for z in courant_core]

L_wm = [2 * x / (o**2) for x, o in zip(Wm, courant_norm)] #calcul de l'inductance avec Wm et I^2 (norme)

R_core = [2 * x / (o**2) for x, o in zip(Pcore, courant_core)] #calcul de R_core avec P_core et I^2

print(f"La valeur de courant_norm est : {courant_norm}")
print(f"La valeur de L_wm est : {L_wm}")
print(f"La valeur de L est : {L}")
print(f"La valeur de R_ohm est : {R_ohm}")
print(f"La valeur de l'impedance Z est : {Impedance}")

R_core_bis = np.zeros(len(R_ohm))
P_coil = np.zeros(len(R_ohm))

for i in range(len(R_ohm)):
    R_core_bis[i] = (R_ohm[i] - R_coil)/(N_turns**2)

print(f"La valeur de R_core est : {R_core}")
print(f"La valeur de R_core_bis est : {R_core_bis}")
for i in range(len(R_ohm)):
    P_coil[i] = 0.5 * R_coil*(np.power(courant_norm[i], 2))

print(f"La valeur de P_coil est : {P_coil}")
print(f"La valeur de Pcoil est : {Pcoil}")

for i in range(len(densite_courant)):
    courant_norm[i] = (2 * Pcore[i]) / R_core[i]
for i in range(len(densite_courant)):
    courant_norm[i] = abs(densite_courant[i])

#print(courant_norm)
#print(R_core)s
#print(abs(courant))
#print(R_coil + R_core)
#print(R_ohm[1])