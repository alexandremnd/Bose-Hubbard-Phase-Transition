import matplotlib
matplotlib.use('TkAgg')

import numpy as np
import matplotlib.pyplot as plt
import textwrap

# Read the data from the file
with open('phase.txt', 'r') as file:
    fixed_param = file.readline().strip()
    data = np.loadtxt(file)

# Extract mu, U, gap_ratio, boson_density, and compressibility
if fixed_param == "J":
    x_label = 'Interaction Strength (U)'
    y_label = 'Chemical Potential (mu)'
    x_values = data[:, 0]
    y_values = data[:, 1]
elif fixed_param == "U":
    x_label = 'Hopping Parameter (J)'
    y_label = 'Chemical Potential (mu)'
    x_values = data[:, 0]
    y_values = data[:, 1]
elif fixed_param == "u":
    x_label = 'Hopping Parameter (J)'
    y_label = 'Interaction Strength (U)'
    x_values = data[:, 0]
    y_values = data[:, 1]
else:
    raise ValueError("Invalid fixed parameter in phase.txt")

gap_ratio = data[:, 2]
boson_density = data[:, 3]
compressibility = data[:, 4]

# Create a grid for x and y
x_unique = np.unique(x_values)
y_unique = np.unique(y_values)
x_grid, y_grid = np.meshgrid(x_unique, y_unique)

# Reshape the data to match the grid
gap_ratio_grid = gap_ratio.reshape(len(y_unique), len(x_unique))
boson_density_grid = boson_density.reshape(len(y_unique), len(x_unique))
compressibility_grid = compressibility.reshape(len(y_unique), len(x_unique))

def wrap_title(title, width=30):
    return "\n".join(textwrap.wrap(title, width))

# Plot the heatmap for gap ratio
plt.figure(figsize=(18, 6))

plt.subplot(1, 3, 1)
plt.contourf(x_grid, y_grid, gap_ratio_grid, levels=50, cmap='viridis')
plt.colorbar(label='Gap Ratio')
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(wrap_title('Gap Ratio with respect to {} and {}'.format(x_label, y_label)), fontsize=10)

# Plot the heatmap for boson density
plt.subplot(1, 3, 2)
plt.contourf(x_grid, y_grid, boson_density_grid, levels=50, cmap='viridis')
plt.colorbar(label='Boson Density')
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(wrap_title('Boson Density with respect to {} and {}'.format(x_label, y_label)), fontsize=10)

# Plot the heatmap for compressibility
plt.subplot(1, 3, 3)
plt.contourf(x_grid, y_grid, compressibility_grid, levels=50, cmap='viridis')
plt.colorbar(label='Compressibility')
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(wrap_title('Compressibility with respect to {} and {}'.format(x_label, y_label)), fontsize=10)

plt.tight_layout()
plt.show()