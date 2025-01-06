import numpy as np
import matplotlib.pyplot as plt

# Data
current = np.array([2.0, 0.9, 0.5, 0.2])
numerical = np.array([1.849, 1.682, 1.54, 1.303829])
theoretical = np.array([1.8388, 1.6768, 1.5337, 1.30335])

magnetic_field = (current * 125) / (10*10**-3)

# Plot
plt.figure(figsize=(8, 5))
plt.scatter(magnetic_field, numerical, color='blue', label=r'Magnetic Flux Density $|\vec{B}_{r_1}|$ with COMSOL')
plt.scatter(magnetic_field, theoretical, color='red', label=r'Saturation Curve $|\vec{B}_{r_1}|$ (Soft Iron)')

# Labels and Legend
plt.xlabel(r'Magnetic field $|\vec{H}_{r_1}|$ (A/m)')
plt.ylabel(r'Magnetic Flux Density $|\vec{B}_{r_1}|$ T)')
plt.title(r'B(H)')
plt.legend()
plt.grid(True)
plt.show()
