import matplotlib
matplotlib.use('TkAgg')

import numpy as np
import matplotlib.pyplot as plt
import textwrap

# Read the data from the file
with open('phase.txt', 'r') as file:
    fixed_param_line = file.readline().strip().split()
    fixed_param = fixed_param_line[0]
    fixed_value = float(fixed_param_line[1])
    data = np.loadtxt(file)

# Extract mu, U, gap_ratio, boson_density, and compressibility
if fixed_param == "J":
    x_label = 'Normalized Interaction Strength (U/J)'
    y_label = 'Normalized Chemical Potential (mu/J)'
    x_values = data[:, 0] / fixed_value
    y_values = data[:, 1] / fixed_value
elif fixed_param == "U":
    x_label = 'Normalized Hopping Parameter (J/U)'
    y_label = 'Normalized Chemical Potential (mu/U)'
    x_values = data[:, 0] / fixed_value
    y_values = data[:, 1] / fixed_value
elif fixed_param == "u":
    x_label = 'Normalized Hopping Parameter (J/mu)'
    y_label = 'Normalized Interaction Strength (U/mu)'
    x_values = data[:, 0] / fixed_value
    y_values = data[:, 1] / fixed_value
else:
    raise ValueError("Invalid fixed parameter in phase.txt")

gap_ratio = data[:, 2]
boson_density = data[:, 3]
compressibility = data[:, 4]

# Create a grid for x and y
ind = np.lexsort((y_values, x_values))
x_values = x_values[ind]
y_values = y_values[ind]
gap_ratio = gap_ratio[ind]

x_unique = np.unique(x_values)
y_unique = np.unique(y_values)
x_grid, y_grid = np.meshgrid(x_unique, y_unique)

# Function to wrap long titles
def wrap_title(title, width=30):
    return "\n".join(textwrap.wrap(title, width))

# Plot the heatmap for gap ratio
plt.figure(figsize=(18, 6))

plt.subplot(1, 3, 1)
contour1 = plt.contourf(x_grid, y_grid, gap_ratio.reshape(len(y_unique), len(x_unique)), levels=50, cmap='viridis')
plt.colorbar(contour1, label='Gap Ratio')
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(wrap_title('Gap Ratio with respect to {} and {}'.format(x_label, y_label)), fontsize=10)

# Plot the heatmap for boson density
plt.subplot(1, 3, 2)
contour2 = plt.contourf(x_grid, y_grid, boson_density.reshape(len(y_unique), len(x_unique)), levels=50, cmap='viridis')
plt.colorbar(contour2, label='Boson Density')
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(wrap_title('Boson Density with respect to {} and {}'.format(x_label, y_label)), fontsize=10)

# Plot the heatmap for compressibility
plt.subplot(1, 3, 3)
contour3 = plt.contourf(x_grid, y_grid, compressibility.reshape(len(y_unique), len(x_unique)), levels=50, cmap='viridis')
plt.colorbar(contour3, label='Compressibility')
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(wrap_title('Compressibility with respect to {} and {}'.format(x_label, y_label)), fontsize=10)

plt.tight_layout()
plt.show()