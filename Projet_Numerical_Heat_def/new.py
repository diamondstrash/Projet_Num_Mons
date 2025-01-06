import matplotlib.pyplot as plt

# Data
current = [2.0, 0.9, 0.5, 0.2]
numerical = [1.8388, 1.6768, 1.5337, 1.3034]
theoretical = [1.849, 1.682, 1.682, 1.3034]

# Plot
plt.figure(figsize=(8, 5))
plt.plot(current, numerical, 'o-', label='Numerical Value', color='blue')
plt.plot(current, theoretical, 's--', label='Theoretical Value', color='red')
plt.xlabel('Current (A)')
plt.ylabel('Magnetic Flux Density (T)')
plt.title('Comparison of Numerical and Theoretical Values')
plt.grid(True)
plt.legend()
plt.show()