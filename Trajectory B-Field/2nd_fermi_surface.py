import numpy as np
import pyvista as pv
from scipy.interpolate import griddata
print(pv.__version__)
# -------------------------------------------------
# 1. LOAD THE DATA
# -------------------------------------------------
# data_file = 'k_surface_fermi_energies_By_0.01_part_1.dat'
# data_file = 'fermi_surface_energies_By_WSM_unfiltered.dat'
data_file = 'arounddirac.dat'
data = np.loadtxt(data_file)
kx = data[:, 0]
ky = data[:, 1]
kz = data[:, 2]
energy = data[:, 4] # Unperturbed band 13
energy13 = data[:, 4] # (Pe)unperturbed

# -------------------------------------------------
# 2. DEFINE FERMI ENERGY
# -------------------------------------------------
unperturbed_min = np.min(energy13) # Minimum of Bottom Conduction Band
# energy_offset = 0.0062  # Energy offset from band minimum (in eV)
energy_offset = 0.001  # Energy offset from band minimum (in eV)

fermi_energy = unperturbed_min + energy_offset
# fermi_energy=4.18903772
print(f"Unperturbed band minimum: {unperturbed_min:.6f} eV")
print(f"Fermi energy set to: {fermi_energy:.6f} eV")

# Filter points near Fermi energy for visualization
# energy_range = np.max(energy) - np.min(energy)
tolerance = 0.0001
lb = fermi_energy - tolerance
ub = fermi_energy + tolerance
mask = (energy >= lb) & (energy <= ub)
kx_filtered = kx[mask]
ky_filtered = ky[mask]
kz_filtered = kz[mask]

# -------------------------------------------------
# 3. CREATE SIMPLE GRID & INTERPOLATE
# -------------------------------------------------
# Keep dimension at 10 as requested
dimension = 150
nx, ny, nz = dimension, dimension, dimension

# # Create a grid with slight expansion 
x_lin = np.linspace(kx_filtered.min(), kx_filtered.max(), nx)
y_lin = np.linspace(ky_filtered.min(), ky_filtered.max(), ny)
z_lin = np.linspace(kz_filtered.min(), kz_filtered.max(), nz)
print(np.min(energy),np.max(energy))
# x_lin = np.linspace(np.min(kx), np.max(kx), nx)
# y_lin = np.linspace(np.min(ky), np.max(ky), ny)
# z_lin = np.linspace(np.min(kz), np.max(kz), nz)
X, Y, Z = np.meshgrid(x_lin, y_lin, z_lin, indexing='ij')
# Basic linear interpolation - use simple approach for speed
print("Performing interpolation...")
points = np.column_stack((kx, ky, kz))
grid_energy = griddata(points, energy, (X, Y, Z), method='linear')
grid_energy = np.nan_to_num(grid_energy, nan=unperturbed_min)


# -------------------------------------------------
# 4. CREATE PYVISTA GRID & EXTRACT ISOSURFACE
# -------------------------------------------------
# Create the grid
grid = pv.StructuredGrid(X, Y, Z)
grid["energy"] = grid_energy.flatten(order='F')

# Extract the basic isosurface - no complicated smoothing
print(f"Extracting iso-surface at energy = {fermi_energy:.6f} eV")
fermi_surface = grid.contour(isosurfaces=[fermi_energy], scalars="energy", method="contour")

# Basic smoothing - minimal iterations to maintain performance
fermi_surface_smooth = fermi_surface.smooth(n_iter=1000, relaxation_factor=0.1,edge_angle=1000)

# -------------------------------------------------
# 5. SIMPLE VISUALIZATION
# -------------------------------------------------
print("Generating visualization...")
# pv.set_plot_theme("dark")
p = pv.Plotter(window_size=[800, 600])

bounds = grid.bounds
background_box = pv.Box(bounds=bounds)

# Add the background box first (so it's behind the data)
p.add_mesh(background_box, color='black', opacity=0.5, show_edges=None, edge_opacity=0.0)
# Add the surface - keep it simple
p.add_mesh(fermi_surface_smooth,
           color='red',
           opacity=0.8,
           specular=0.8,
           specular_power=15,
           smooth_shading=True,
           show_scalar_bar=False)
# Custom axis range
# show_bounds returns a vtkCubeAxesActor object
axes_actor = p.show_bounds(
    grid='back',
    location='outer',   # draw the grid behind the data
    color='Black',
    show_xaxis=True,
    show_yaxis=True,
    show_zaxis=True,  # white ticks/lines on black background
    ticks='outside',   # automatically choose tick positions
    all_edges=False, # show ticks on all edges
    fmt='%.2f',    # Format tick labels to 2 decimal places
    xtitle="k_x", 
    ytitle="k_y",
    ztitle="k_z",
    font_size=12,
    font_family="arial",
    n_zlabels=3,
)


# Add a bounding box to show the k-space domain
outline = grid.outline()
p.add_mesh(outline, color='gray', opacity=0.25, line_width=2)
# Simple coordinate axes
p.add_light(pv.Light(position=(1, 1, 1), focal_point=(0, 0, 0), 
                    color=[1, 1, 1], intensity=0.8))
p.add_light(pv.Light(position=(-1, -1, -1), focal_point=(0, 0, 0), 
                    color=[0.5, 0.5, 0.7], intensity=0.3))
# Add title
# p.add_title(f"Fermi Surface at E = {fermi_energy:.4f} eV")
# p.window_size = [2000, 1500]  # Increase resolution


# Show the plot
p.show()