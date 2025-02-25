import numpy as np
from numpy import sin, cos, pi
# from skimage import measure
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
from scipy.interpolate import griddata

# -------------------------------------------------
# 1. LOAD THE DATA
# -------------------------------------------------
# Replace 'data.txt' with the path to your file
# If your file has 5 columns, adjust slicing accordingly (e.g., kx, ky, kz, E, Fermi).
data = np.loadtxt('k_surface_fermi_energies.dat')
kx = data[:, 0]
ky = data[:, 1]
kz = data[:, 2]
energy = data[:, 3]
energy13 = data[:,4]
# If your file provides a per-row Fermi energy in the 4th column and the actual energy in the 5th,
# just swap the indexing. The key idea is to extract each column correctly.

# -------------------------------------------------
# 2. DEFINE FERMI ENERGY AND FILTER
# -------------------------------------------------
# fermi_energy = 4.18903772  # Provided Fermi energy from wanr2k.f90 file

lowest_energy = np.min(energy13)
#fermi_energy = 4.198475999999999  # Provided Fermi energy set from lowest_energy + 0.02
fermi_energy = lowest_energy
lower_bound = fermi_energy - 0.012
upper_bound = fermi_energy + 0.012
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
nx, ny, nz = 10, 10, 10

# Create a regular grid covering the region of filtered k-space.
x_lin = np.linspace(kx_filtered.min(), kx_filtered.max(), nx)
y_lin = np.linspace(ky_filtered.min(), ky_filtered.max(), ny)
z_lin = np.linspace(kz_filtered.min(), kz_filtered.max(), nz)
X, Y, Z = np.meshgrid(x_lin, y_lin, z_lin, indexing='ij')

# Interpolate the energy values onto the grid.
points = np.column_stack((kx, ky, kz))
grid_energy = griddata(points, energy, (X, Y, Z), method='linear')

# If there are NaNs (areas not covered by data), you might fill them.
grid_energy = np.nan_to_num(grid_energy, nan=lowest_energy)

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
contours = grid.contour(isosurfaces=[fermi_energy], scalars="energy")
# inner_surface = grid.contour(isosurfaces=[lower_bound], scalars="energy")

# -------------------------------------------------
# 6. VISUALIZE THE RESULT
# -------------------------------------------------
pv.set_plot_theme("document")
p = pv.Plotter()
# Color the mesh by its z-coordinate (or any scalar you prefer)
# p.add_mesh(inner_surface, opacity=0.5, scalars=inner_surface.points[:, 2], show_scalar_bar=True)
p.add_mesh(contours, color="red", show_scalar_bar=True, opacity=0.4, scalars=contours.points[:, 2])#, show_scalar_bar=True)
p.add_axes()
p.add_title("Iso-Surface at Fermi Energy")
p.show()