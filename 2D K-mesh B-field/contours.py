import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.widgets import Slider
# Load the data
data = np.loadtxt('Energy_part_0.5B_np500_1.dat')

kx = data[:, 0]
ky = data[:, 1]
energy = data[:, 2]
energy -= np.min(energy)

trian = tri.Triangulation(kx,ky)
#energies = [0,  0.01, 0.02, 0.04, 0.170, 0.190]
# Plot
fig, ax = plt.subplots(figsize=(6,6))
plt.subplots_adjust(bottom=0.25)  # Leave space for slider

ax.set_xlim([-0.1, 0.1])      # Set kx range
ax.set_ylim([-0.1, 0.1])      # Set ky range

# Plot the contour
initial_level = np.min(energy) + 0.1 * (np.max(energy) - np.min(energy))
contour = ax.tricontour(trian, energy, levels=[initial_level], linewidths=1.5, cmap="viridis")
lines = contour
plt.xlabel('$k_x$ ($A^{-1}$)')
plt.ylabel('$k_y$ ($A^{-1}$)')
plt.title('Contour Plot of Energy')
#plt.colorbar(contour, label='Energy (eV)')
# for i, energy in enumerate(energies.T):
#     ax.plot_trisurf(kx, ky, energy, cmap='viridis', edgecolor='none', alpha=0.7, label=f'Column {i+4}')
# ax.view_init(elev=30, azim=135)
ax_slider = plt.axes([0.2, 0.1, 0.6, 0.03])  # Position of slider
energy_slider = Slider(ax_slider, 'Energy Level', np.min(energy), 0.1, valinit=initial_level)
# Update function
def update(val):
    # Remove all previous contours
    while ax.collections:
        ax.collections[-1].remove()

    # Get new energy level from slider
    selected_level = energy_slider.val

    # Draw new contour at selected energy level
    ax.tricontour(trian, energy, levels=[selected_level], linewidths=2, colors='blue')

    # Refresh the figure
    plt.draw()

# plt.colorbar(contour, label='Energy (eV)')
energy_slider.on_changed(update)
plt.show()
