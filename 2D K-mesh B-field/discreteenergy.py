import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy.interpolate import griddata

# Load the data from file
# data = np.loadtxt('NewDataWSMUnperturbed1.dat')
data = np.loadtxt('NewDataWSMPerturbed1.dat')

# Extract columns
kx = data[:, 0]
ky = data[:, 1]
energy = data[:, 10]
energy -= np.min(energy)
# 0.001 ev = 1 meV, so 10meV = 0.010, 170 meV = 0.170 eV
# Create a regular grid for the contour plot
grid_x, grid_y = np.mgrid[min(kx):max(kx):1000j, min(ky):max(ky):1000j]

# Interpolate the irregular data onto a regular grid
grid_z = griddata((kx, ky), energy, (grid_x, grid_y))

# Define custom energy levels for the contour plot
# You can modify these values to suit your specific energy ranges
# energy_levels = np.linspace(min(energy), max(energy), 10)  # 10 evenly spaced levels
# Alternatively, specify exact levels you want:
energy_levels = [0.001, 0.01, 0.0226, 0.068, 0.15, 0.2215]
# energy_levels = [0.001, 0.01, 0.0226, 0.05, 0.12, 0.2]

custom_cmap = ['black', 'purple', 'blue', 'orange', 'green', 'red']

# Create figure and axes
fig, ax = plt.subplots(figsize=(10, 8))
ax.set_xlim([-0.1, 0.1])
ax.set_ylim([-0.1, 0.1])
# Simply set fewer ticks with set_xticks and set_yticks
ax.set_xticks([-0.1, -0.05, 0, 0.05, 0.1])
ax.set_yticks([-0.1, -0.05, 0, 0.05, 0.1])
ax.tick_params(axis='both', which='major', labelsize=16, pad=10)  # Increased font size for tick labels

# Create the contour plot with custom levels
contour = ax.contour(grid_x, grid_y, grid_z, levels=energy_levels, colors=custom_cmap, linewidths = 4)

# energy_labels = ['0', '10', '20', '50', '120', '200'] #Unperturbed
energy_labels = ['0', '10', '20', '70', '150', '220'] #Perturbed

cbar = fig.colorbar(contour, ax=ax, location='top', orientation='horizontal', shrink=0.6)

# Set tick labels and increase their font size
cbar.set_ticklabels(energy_labels)
cbar.ax.tick_params(labelsize=20)  # Adjust font size as needed

# Make the colorbar lines thicker and more visible
# For contour colorbars, set the tick lines and extend the colorbar height
# Make the colorbar lines thicker and set their colors to match the contour lines
# cbar.ax.tick_params(width=3, length=15)  # Thicker, longer ticks
# for tickline, color in zip(cbar.ax.xaxis.get_ticklines(), custom_cmap * 20):  # *2 in case of both sides
#   tickline.set_linewidth(3)
#   tickline.set_color(color)
# cbar.outline.set_linewidth(3)  # Thicker colorbar outline
cbar.set_label(r'$E-E_{\mathrm{CBM}}$ (meV)', fontsize=25, labelpad=10)

# Add labels and title
ax.set_xlabel(r'$k_x$  (Å$^{-1}$)', fontsize=25)
ax.set_ylabel(r'$k_y$  (Å$^{-1}$)',  fontsize=25)
# ax.set_title('Energy Contour Plot')
ax.set_xticks
# Make the plot more aesthetically pleasing
# plt.grid(True, linestyle='--', alpha=1)
ax.set_aspect('equal')  # Equal aspect ratio
# Save the plot in high resolution before displaying
plt.tight_layout()
plt.savefig('energy_contourPERT.png', dpi=1000, bbox_inches='tight')  # You can change filename and dpi as needed
plt.show()