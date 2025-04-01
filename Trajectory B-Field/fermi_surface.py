import numpy as np
import pyvista as pv
from scipy.interpolate import griddata

# -------------------------------------------------
# 1. LOAD THE DATA
# -------------------------------------------------
# Replace 'data.txt' with the path to your file
# If your file has 5 columns, adjust slicing accordingly (e.g., kx, ky, kz, E, Fermi).
data_file = 'k_surface_fermi_energies_By_WSM_part_1.dat'
data = np.loadtxt(data_file)
# band = 4
kx = data[:, 0]
ky = data[:, 1]
kz = data[:, 2]
energy = data[:, 3]
energy13 = data[:, 4]
# If your file provides a per-row Fermi energy in the 4th column and the actual energy in the 5th,
# just swap the indexing. The key idea is to extract each column correctly.

# -------------------------------------------------
# 2. DEFINE FERMI ENERGY WITH PHYSICAL MEANING
# -------------------------------------------------
# Set Fermi energy based on unperturbed band minimum (CBM)
unperturbed_min = np.min(energy13)
energy_offset = 0.02  # Energy offset from band minimum (in eV)
fermi_energy = unperturbed_min + energy_offset

print(f"Unperturbed band minimum: {unperturbed_min:.6f} eV")
print(f"Fermi energy set to: {fermi_energy:.6f} eV (+{energy_offset} eV from band minimum)")

energy_range = np.max(energy) - np.min(energy)
tolerance = energy_range * 0.001  # Adaptive tolerance (0.1% of total range)
lower_bound = fermi_energy - tolerance
upper_bound = fermi_energy + tolerance
mask = (energy >= lower_bound) & (energy <= upper_bound)

# Filter the scattered data (just for a quick 3D scatter)
kx_filtered = kx[mask]
ky_filtered = ky[mask]
kz_filtered = kz[mask]
energy_filtered = energy[mask]
# -------------------------------------------------
# 3. INTERPOLATE THE FILTERED DATA ONTO A 3D GRID
# -------------------------------------------------
# Define the grid resolution. Increase nx, ny, nz for finer detail.
dimension = 10
nx, ny, nz = dimension, dimension, dimension

# # Create a regular grid covering the region of filtered k-space.
# x_lin = np.linspace(kx.min(), kx.max(), nx)
# y_lin = np.linspace(ky.min(), ky.max(), ny)
# z_lin = np.linspace(kz.min(), kz.max(), nz)

expand_factor = 0.05
x_min, x_max = kx.min() - expand_factor*abs(kx.min()), kx.max() + expand_factor*abs(kx.max())
y_min, y_max = ky.min() - expand_factor*abs(ky.min()), ky.max() + expand_factor*abs(ky.max())
z_min, z_max = kz.min() - expand_factor*abs(kz.min()), kz.max() + expand_factor*abs(kz.max())

x_lin = np.linspace(x_min, x_max, nx)
y_lin = np.linspace(y_min, y_max, ny)
z_lin = np.linspace(z_min, z_max, nz)
X, Y, Z = np.meshgrid(x_lin, y_lin, z_lin, indexing='ij')

# Interpolate the energy values onto the grid.
print("Performing high-quality cubic interpolation (this may take a moment)...")
points = np.column_stack((kx, ky, kz))
grid_energy = griddata(points, energy, (X, Y, Z), method='linear')

# If there are NaNs (areas not covered by data), you might fill them.
if np.any(np.isnan(grid_energy)):
    print(f"Filling {np.sum(np.isnan(grid_energy))} NaN values with nearest-neighbor interpolation")
    grid_energy_nearest = griddata(points, energy, (X, Y, Z), method='nearest')
    grid_energy = np.where(np.isnan(grid_energy), grid_energy_nearest, grid_energy)

# -------------------------------------------------
# 4. CREATE A STRUCTURED GRID WITH PYVISTA
# -------------------------------------------------
# PyVista expects the coordinate arrays to have the same shape.
grid = pv.StructuredGrid(X, Y, Z)
# Attach the interpolated energy field as point data.
# Note: Flatten in Fortran order ('F') to match PyVista's internal ordering.
grid["energy"] = grid_energy.flatten(order='F')

# -------------------------------------------------
# 5. EXTRACT THE ISO-SURFACE (CONTOUR)
# -------------------------------------------------
# For example, we extract the isosurface at the Fermi energy.
print(f"Extracting iso-surface at energy = {fermi_energy:.6f} eV")
fermi_surface = grid.contour(isosurfaces=[fermi_energy], scalars="energy")
# inner_surface = grid.contour(isosurfaces=[lower_bound], scalars="energy")
fermi_surface_smooth = fermi_surface.smooth(n_iter=100, relaxation_factor=0.2)
# -------------------------------------------------
# 6. VISUALIZE THE RESULT
# -------------------------------------------------
print("Generating high-quality visualization...")
pv.set_plot_theme("document")
p = pv.Plotter(window_size=[1200, 1000])

# Add the original points as small dots to verify accuracy
p.add_points(np.column_stack((kx_filtered, ky_filtered, kz_filtered)),
             color='black', point_size=3, opacity=0.3)
p.set_background("black")
# Add the Fermi surface with enhanced rendering
p.add_mesh(fermi_surface_smooth, 
           color='red',                # Try 'blue', 'green', or other colors
           opacity=0.7,               # Adjust for best visibility
           smooth_shading=True,       # Enable Phong shading for smoothness
           specular=0.5,              # Add specular highlight for 3D effect
           specular_power=15,         # Adjust specular power
           show_scalar_bar=False)     # Toggle to True if you want a color bar
p.show_bounds(grid='back',               # draw the grid behind the data
                color='white',             # white lines/ticks against black background
                outline=True,              # show the bounding box outline
                ticks='auto',              # automatically choose tick locations
                all_edges=True,            # show all edges with ticks
                font_size=12,              # increase label font size
                linewidth=2,               # slightly thicker lines
                xtitle='kₓ (Å⁻¹)',
                ytitle='kᵧ (Å⁻¹)',
                ztitle='k𝓏 (Å⁻¹)')
# Add coordinate axes with proper labels
p.add_axes(xlabel='k_x', ylabel='k_y', zlabel='k_z', 
           line_width=4, labels_off=False)

# Add a bounding box to show the k-space domain

outline = grid.outline()
p.add_mesh(outline, color='gray', opacity=0.8, line_width=3)


# Add proper lighting for better 3D perception
p.add_light(pv.Light(position=(1, 1, 1), focal_point=(0, 0, 0), 
                    color=[1, 1, 1], intensity=0.8))
p.add_light(pv.Light(position=(-1, -1, -1), focal_point=(0, 0, 0), 
                    color=[0.5, 0.5, 0.7], intensity=0.3))

# Add title and information
p.add_title(f"Fermi Surface at E = {fermi_energy:.4f} eV", font_size=16)
p.show()