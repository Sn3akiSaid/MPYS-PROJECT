import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcParams
import matplotlib as mpl

# Create custom colormap (equivalent to gnuplot's palette)
colors = [(1, 0, 0), (1, 1, 1), (0, 0, 1)]  # red, white, blue
cmap_name = 'bwr'
# cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

# Create figure with the specified dimensions
# Convert inches to inches (keeping the 5x4.5 size)
fig, ax = plt.subplots(figsize=(5, 4.5), dpi=100)

# Load data
data = np.loadtxt('band_partition_LAH_WSM_1.dat')

# Extract columns based on gnuplot's usage
# x1 = data[::5, 0]  # First column (k_x)
# y1 = data[::5, 1]  # Second column (E-E_F for first plot)
# c1 = data[::5, 2]  # Fourth column (color data for first plot)

# x2 = data[:, 0]  # First column (k_x again)
# y2 = data[:, 3]  # Sixth column (E-E_F for second plot)
# c2 = data[:, 4]  # Eighth column (color data for second plot)
x1 = data[::5, 0]  # First column (k_x)
y1 = data[::5, 1]  # Second column (E-E_F for first plot)
c1 = data[::5, 3]  # Fourth column (color data for first plot)

x2 = data[:, 0]  # First column (k_x again)
y2 = data[:, 5]  # Sixth column (E-E_F for second plot)
c2 = data[:, 7]  # Eighth column (color data for second plot)
# Normalize color values to match gnuplot's range
# norm = mpl.colors.Normalize(vmin=-1, vmax=1)

ax.scatter(x1, y1, c=c1, cmap=cmap_name, linestyle='--', s=6)
ax.scatter(x2, y2, c=c2, cmap=cmap_name, s=6)
# Plot the second dataset with solid lines
# for i in range(len(x2)-1):
#     # Use color from the colormap based on c2 value
#     color = cm(norm(c2[i]))
#     ax.plot(x2[i:i+2], y2[i:i+2], color=color, s=2)

# Set axis ranges to match gnuplot
ax.set_xlim(-0.14, 0.14)
ax.set_ylim(-0.5, 0.5)

# Set axis labels with subscripts and superscripts to match gnuplot
ax.set_xlabel(r'$k_{x}$ (Å$^{-1}$)')
ax.set_ylabel(r'$E-E_{F}$ (eV)')

# Adjust the position of the y-label to match the offset in gnuplot
# ax.yaxis.set_label_coords(-0.12, 0.5)  # Approximating the offset 1.2,0 from gnuplot

# Add a colorbar to show the color scale
# sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
# sm.set_array([])
# cbar = plt.colorbar(sm, ax=ax)
# cbar.set_ticks([-1, 0, 1])

# Adjust the layout and save
plt.tight_layout()
# plt.savefig('DatayLAHWSM.pdf', format='pdf', bbox_inches='tight')

# Show the plot (optional, comment out for headless environments)
plt.show()