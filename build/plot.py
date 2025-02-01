import matplotlib
matplotlib.use('TkAgg')

import numpy as np
import matplotlib.pyplot as plt

# Read the data from the file
data = np.loadtxt('phase.txt')

# Extract mu, U, gap_ratio, boson_density, and compressibility
mu_values = data[:, 0]
U_values = data[:, 1]
gap_ratio = data[:, 2]
boson_density = data[:, 3]
compressibility = data[:, 4]

# Create a grid for mu and U
mu_unique = np.unique(mu_values)
U_unique = np.unique(U_values)
mu_grid, U_grid = np.meshgrid(mu_unique, U_unique)

# Reshape the data to match the grid
gap_ratio_grid = gap_ratio.reshape(len(U_unique), len(mu_unique))
boson_density_grid = boson_density.reshape(len(U_unique), len(mu_unique))
compressibility_grid = compressibility.reshape(len(U_unique), len(mu_unique))

# Plot the heatmap for gap ratio
plt.figure(figsize=(18, 6))

plt.subplot(1, 3, 1)
plt.contourf(mu_grid, U_grid, gap_ratio_grid, levels=50, cmap='viridis')
plt.colorbar(label='Gap Ratio')
plt.xlabel('Chemical Potential (mu)')
plt.ylabel('Interaction Strength (U)')
plt.title('Gap Ratio with respect to U and mu')

# Plot the heatmap for boson density
plt.subplot(1, 3, 2)
plt.contourf(mu_grid, U_grid, boson_density_grid, levels=50, cmap='viridis')
plt.colorbar(label='Boson Density')
plt.xlabel('Chemical Potential (mu)')
plt.ylabel('Interaction Strength (U)')
plt.title('Boson Density with respect to U and mu')

# Plot the heatmap for compressibility
plt.subplot(1, 3, 3)
plt.contourf(mu_grid, U_grid, compressibility_grid, levels=50, cmap='viridis')
plt.colorbar(label='Compressibility')
plt.xlabel('Chemical Potential (mu)')
plt.ylabel('Interaction Strength (U)')
plt.title('Compressibility with respect to U and mu')

plt.tight_layout()
plt.show(block=True)
exit()