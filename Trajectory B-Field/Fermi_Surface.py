
import vtk
import numpy as np
import pyvista as pv
from scipy.interpolate import griddata
print(pv.__version__)
print(pv.Report())
# pv.global_theme.return_cpos = True
# -------------------------------------------------
# 1. LOAD THE DATA
# -------------------------------------------------
# data_file = 'k_surface_fermi_energies_By_0.01_part_1.dat'

# data_file = 'k_surface_fermi_energies.dat'


# data_file = 'berryinWSM_11.dat'
# data_file = 'berryinWSM_perturbed16.dat'
# data_file = "NEWberryinWSM_perturbedBigger1.dat" #GOOD PHASE
# data_file = "NEWberryinWSM_perturbedBigger0749.dat"
# data_file = "BerryPertTVB_15.dat"
data_file = "BerryUnpertTVB_5.dat"

data = np.loadtxt(data_file)
kx = data[:, 0]
ky = data[:, 1]
kz = data[:, 2]
omegax = data [:,3]
omegay = data [:,4]
omegaz = data [:,5]
absomega = data [:,6]
energy = data[:, 7]
energy13 = data[:, 7]

#FERMI
# kx = data[:, 0]
# ky = data[:, 1]
# kz = data[:, 2]
# energy = data[:, 4] # Unperturbed band 13
# energy13 = data[:, 4] # (Pe)unperturbed

# -------------------------------------------------
# 2. DEFINE FERMI ENERGY WITH PHYSICAL MEANING
# -------------------------------------------------
unperturbed_min = np.min(energy13)#5.577098
# energy_offset = 0.018
# energy_offset = 0.006 #pert # Energy offset from band minimum (in eV)
energy_offset = 0.01
# energy_offset = 0.00624
fermi_energy = unperturbed_min + energy_offset
# fermi_energy=4.18903772
print(f"Unperturbed band minimum: {unperturbed_min:.6f} eV")
print(f"Fermi energy set to: {fermi_energy:.6f} eV")

# Filter points near Fermi energy for visualization
# energy_range = np.max(energy) - np.min(energy)
tolerance = 0.5
lb = fermi_energy - tolerance
ub = fermi_energy + tolerance
mask = (energy >= lb) & (energy <= ub)
kx_filtered = kx[mask]
ky_filtered = ky[mask]
kz_filtered = kz[mask]
omegax_filtered = omegax[mask]# / absomega[mask] 
omegay_filtered = omegay[mask]# / absomega[mask]
omegaz_filtered = omegaz[mask]# / absomega[mask]
# -------------------------------------------------
# 3. CREATE SIMPLE GRID & INTERPOLATE
# -------------------------------------------------
# Keep dimension at 10 as requested
dimension = 100
nx, ny, nz = dimension, dimension, dimension

# # Create a grid with slight expansion
print(kx_filtered.min(), kx_filtered.max())
# x_lin = np.linspace(kx_filtered.min(), kx_filtered.max(), nx)
# y_lin = np.linspace(ky_filtered.min(), ky_filtered.max(), ny)
# z_lin = np.linspace(kz_filtered.min(), kz_filtered.max(), nz)

x_lin = np.linspace(-0.08, 0.08, nx)
y_lin = np.linspace(-0.08, 0.08,  ny)
z_lin = np.linspace(0.42,0.51, nz)
X, Y, Z = np.meshgrid(x_lin, y_lin, z_lin, indexing='ij')
# Basic linear interpolation - use simple approach for speed
print("Performing interpolation...")
points = np.column_stack((kx, ky, kz))
grid_energy = griddata(points, energy, (X, Y, Z), method='nearest')
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
fermi_surface_smooth = fermi_surface.smooth(n_iter=10000, feature_angle=10000, relaxation_factor=0.01,edge_angle=10000)

# -------------------------------------------------
# 5. SIMPLE VISUALIZATION
# -------------------------------------------------
print("Generating visualization...")
# pv.set_plot_theme("dark")
p = pv.Plotter(window_size=[1440,1080], 
               line_smoothing=True,
               polygon_smoothing=True, )
            #    border=True,
            #    border_color='red')

bounds = grid.bounds
background_box = pv.Box(bounds=bounds)

# Add the background box first (so it's behind the data)
p.add_mesh(background_box, color='black', opacity=0.2, show_edges=None, edge_opacity=0.0)
# Add the surface - keep it simple
p.add_mesh(fermi_surface_smooth,
           color='red',
           opacity=1,
           specular=0.9,
           specular_power=50,
           smooth_shading=True,
           show_scalar_bar=False,
           silhouette=True,
           lighting=True,
        #    ambient=0.4,
           diffuse=0.9,
           )

# Custom axis range
# show_bounds returns a vtkCubeAxesActor object
axes_actor = p.show_bounds(
    grid='back',
    # axes_ranges=(-0.06, 0.06, -0.06, 0.06, 0.42, 0.49),
    # bounds=[-0.06, 0.06, -0.06, 0.06, 0.42, 0.49],
    location='outer',   # draw the grid behind the data
    color='Black',
    show_xaxis=True,
    show_yaxis=True,
    show_zaxis=True,  # white ticks/lines on black background
    ticks='outside',   # automatically choose tick positions
    all_edges=False, # show ticks on all edges
    # corner_factor=1,
    fmt='%.2f',    # Format tick labels to 2 decimal places
    # xtitle=r"$k_x(\AA^{-1})$", 
    # ytitle=r"$k_y(\AA^{-1})$",
    # ztitle=r"$k_z(\AA^{-1})$",
    xtitle="", 
    ytitle="",
    ztitle="",
    
    font_size=20,
    font_family="arial",
    n_xlabels=3,
    n_ylabels=3,
    n_zlabels=3,
    # padding=0.01
)
# After creating your bounds with show_bounds()
# Get the cube axes actor
cube_axes_actor = p.renderer.cube_axes_actor
cube_axes_actor.label_offset = 40.0
# cube_axes_actor.SetUseTextActor3D(True)
# cube_axes_actor.GetTitleTextProperty(0).SetFontSize(1)  # X axis title
# cube_axes_actor.GetTitleTextProperty(1).SetFontSize(1)  # Y axis title
# cube_axes_actor.GetTitleTextProperty(2).SetFontSize(1)  # Z axis title
# Increase the tick offset
# This pushes the tick labels further away from the axes
# cube_axes_actor.SetXLabelOffset(20)  # Default is around 2
# cube_axes_actor.SetYLabelOffset(20)
# cube_axes_actor.SetZLabelOffset(20)

# # If needed, also adjust the label positions
# cube_axes_actor.SetTickLocationToOutside()
# p.reset_camera()
# p.render()
# bounds = grid.bounds
# x_center = (bounds[0] + bounds[1])/2
# y_center = (bounds[2] + bounds[3])/2
# z_center = (bounds[4] + bounds[5])/2

# # Add 3D text for axis labels
# p.add_point_labels(
#     [[x_center, bounds[3]+bounds[3], bounds[4]], 
#      [bounds[1]+bounds[1], y_center, bounds[4]],
#      [bounds[1]+bounds[1], bounds[2], z_center]],
#     ["$k_x$", "$k_y$", "$k_z$"],
#     font_size=12,
#     always_visible=True,
#     shadow=False,
#     shape=None,
#     point_size=0,  # Hide the points
# )
# p.add_text("$k_x$", position='lower_edge', font_size=8, color='black')
# p.add_text("$k_y$", position='left_edge', font_size=8, color='black')
# p.add_text("$k_z$", position='upper_left', font_size=8, color='black')
# FONT SIZES
# axes_actor.GetTitleTextProperty(0).SetFontSize(10)  # X-axis title
# axes_actor.GetTitleTextProperty(1).SetFontSize(10)  # Y-axis title
# axes_actor.GetTitleTextProperty(2).SetFontSize(10)  # Z-axis title

# # TICKS
# axes_actor.GetLabelTextProperty(0).SetFontSize(10)  # X-axis tick labels
# axes_actor.GetLabelTextProperty(1).SetFontSize(10)  # Y-axis tick labels
# axes_actor.GetLabelTextProperty(2).SetFontSize(10)  # Z-axis tick labels

# Add a bounding box to show the k-space domain
# outline = grid.outline()
# p.add_mesh(outline, color='gray', opacity=0.25, line_width=0)
# Simple coordinate axes
p.add_light(pv.Light(position=(1, 1, 1), focal_point=(0, 0, 0), 
                    color=[1, 1, 1], intensity=0.8))
p.add_light(pv.Light(position=(-1, -1, -1), focal_point=(0, 0, 0), 
                    color=[0.5, 0.5, 0.7], intensity=0.3))
# Add title
# p.add_title(f"Fermi Surface at E = {fermi_energy:.4f} eV")
# p.window_size = [2000, 1500]  # Increase resolution
# -------------------------------------------------
# 6. ADD OMEGA VECTORS TO VISUALIZATION
# -------------------------------------------------
print("Adding omega vectors to visualization...")

# We need to identify points that are actually on the Fermi surface
# Get the first few layers of surface points from the fermi_surface_smooth
surface_points = np.array(fermi_surface_smooth.points)

# Create a KDTree for efficient nearest-neighbor lookup
from scipy.spatial import KDTree
tree = KDTree(surface_points)

# For each filtered point, check if it's close to the Fermi surface
# Only keep points that are within a small distance of the surface
filtered_points = np.column_stack((kx_filtered, ky_filtered, kz_filtered))
# max_distance = 0.005  # Maximum distance to be considered "on" the surface - adjust as needed
max_distance = 0.012
# Query the KDTree to find distances to nearest surface points
distances, _ = tree.query(filtered_points, k=1)
surface_mask = distances <= max_distance

# Use the mask to select only points near the surface
surface_kx = kx_filtered[surface_mask]
surface_ky = ky_filtered[surface_mask]
surface_kz = kz_filtered[surface_mask]
surface_omegax = omegax_filtered[surface_mask]
surface_omegay = omegay_filtered[surface_mask]
surface_omegaz = omegaz_filtered[surface_mask]

print(f"Original filtered points: {len(kx_filtered)}")
print(f"Points on Fermi surface: {len(surface_kx)}")

# Sample points if there are too many
# Create a uniform 3D grid for sampling
def grid_sample_3d_with_ranges(x, y, z, vx, vy, vz, num_bins, kz_ranges=None):
    """
    Sample points evenly throughout 3D space using grid cells
    with optional filtering for specific kz ranges
    
    Parameters:
    -----------
    kz_ranges : list of tuples
        List of (min, max) ranges for kz to include
        e.g. [(0.42, 0.44), (0.46, 0.48), (0.50, 0.51)]
        If None, all kz values are included
    """
    # First filter by kz ranges if specified
    if kz_ranges is not None:
        # Create a mask for all points in any of the kz ranges
        kz_mask = np.zeros_like(z, dtype=bool)
        for kz_min, kz_max in kz_ranges:
            kz_mask = kz_mask | ((z >= kz_min) & (z <= kz_max))
        
        # Filter all arrays by this mask
        if not np.any(kz_mask):
            print("No points found in the specified kz ranges")
            return (np.array([]), np.array([]), np.array([]),
                   np.array([]), np.array([]), np.array([]))
        
        x = x[kz_mask]
        y = y[kz_mask]
        z = z[kz_mask]
        vx = vx[kz_mask]
        vy = vy[kz_mask]
        vz = vz[kz_mask]
    
    # Continue with normal grid sampling on filtered points
    # Create 3D histogram bins
    x_bins = np.linspace(np.min(x), np.max(x), num_bins)
    y_bins = np.linspace(np.min(y), np.max(y), num_bins)
    z_bins = np.linspace(np.min(z), np.max(z), num_bins)
    
    # Lists to store selected points
    selected_x, selected_y, selected_z = [], [], []
    selected_vx, selected_vy, selected_vz = [], [], []
    
    # For each 3D grid cell, select one point (if any exist)
    for i in range(len(x_bins)-1):
        for j in range(len(y_bins)-1):
            for k in range(len(z_bins)-1):
                # Find points in this cell
                mask = (
                    (x >= x_bins[i]) & (x < x_bins[i+1]) &
                    (y >= y_bins[j]) & (y < y_bins[j+1]) &
                    (z >= z_bins[k]) & (z < z_bins[k+1])
                )
                indices = np.where(mask)[0]
                
                if len(indices) > 0:
                    # Take the first point in this cell
                    idx = indices[0]
                    selected_x.append(x[idx])
                    selected_y.append(y[idx])
                    selected_z.append(z[idx])
                    selected_vx.append(vx[idx])
                    selected_vy.append(vy[idx])
                    selected_vz.append(vz[idx])
    
    return (np.array(selected_x), np.array(selected_y), np.array(selected_z),
            np.array(selected_vx), np.array(selected_vy), np.array(selected_vz))

# Adjust num_bins to control density (higher = more vectors)
# For 5000 vectors, try num_bins around 20-30
# kz_ranges = [
#     (0.42, 0.45),  # First range
#     (0.455, 0.465),  # Second range
#     (0.465, 0.51)   # Third range
# ]
num_bins = 16#pert
# num_bins = 22
s_kx, s_ky, s_kz, s_ox, s_oy, s_oz = grid_sample_3d_with_ranges(
    surface_kx, surface_ky, surface_kz,
    surface_omegax, surface_omegay, surface_omegaz,
    num_bins,
    # kz_ranges=kz_ranges
)

# Continue with your existing code using these sampled points
surface_points = np.column_stack((s_kx, s_ky, s_kz))
surface_point_cloud = pv.PolyData(surface_points)
surface_vectors = np.column_stack((s_ox, s_oy, s_oz))
surface_point_cloud["omega_vectors"] = surface_vectors

# Create the glyphs
# Normalize vectors to have uniform length, but keep direction
vector_norms = np.sqrt(np.sum(surface_vectors**2, axis=1))
scale_factor = 0.0035  # Adjust as needed

custom_arrow = pv.Arrow(
    start=(-0.5, 0, 0),         # Starting point of the arrow
    direction=(1, 0, 0),      # Direction of the arrow
    tip_length=0.65,          # Length of the tip as a fraction of the total length
    tip_radius=0.19,           # Radius of the tip
    tip_resolution=50,        # Number of faces around the tip
    shaft_radius=0.1,        # Radius of the shaft
    shaft_resolution=50,      # Number of faces around the shaft
    scale=2.5                # Overall scaling factor
)

# Create the glyphs
glyphs = surface_point_cloud.glyph(
    orient="omega_vectors",
    scale=False,
    factor=scale_factor,
    geom=custom_arrow,
    tolerance=0.0,
    # color_mode='vector'
)

# Add the glyphs to the plot
p.add_mesh(glyphs, 
    color=[0,192,255],
    # color="#4ebcff",
    show_scalar_bar=False,
    name="omega_vectors",
    specular=0.8,
    specular_power=50,
    smooth_shading=True,
      # Add these parameters:
    lighting=True,                 # Enhance 3D appearance with lighting
    opacity=1,                   # Slightly transparent for better visibility
    line_width=3.0,                # Thicker lines when viewing wireframes
    render_points_as_spheres=True, # Better quality point rendering
    ambient=0.4,                   # Increase ambient light to see vectors better
    diffuse=0.8,                   # Good diffuse scattering for 3D appearance
    pickable=False,                 # Prevent accidental selection
    # show_edges=True,
    # edge_color="black",
    silhouette={'color': 'black',       # Silhouette color
                'line_width': 2.8,      # Width of the silhouette lines
                'opacity': 1,         # Transparency of silhouette (0-1)
                'feature_angle': 45.0,  # Display edges exceeding this angle
                'decimate': 0.5         # Decimation level (reduces mesh complexity)
                }
)

# # Optional: Add a legend
# # p.add_legend([("Fermi Surface", "red"), ("Berry Curvature", "blue")])

# # Show the plot
p.enable_anti_aliasing('ssaa')
# p.camera.azimuth=135 #perturbed wsm
p.camera.azimuth=170
# p.camera.zenith=50
current_focal_point = p.camera.focal_point
new_focal_point = [
    current_focal_point[0],# - 0.01,  # Shift left (-x direction)
    current_focal_point[1],# - 0.001,  # Shift up (+y direction)
    current_focal_point[2] - 0.02         # Keep same z
]
p.camera.focal_point = new_focal_point
p.show(return_cpos=True,return_img=True)
# p.screenshot('3dFermiplotTRIVUNPert.png', transparent_background=True)

# p.screenshot(
# filename=None,
# transparent_background=None,
# return_img=True,
# window_size=None,
# scale=None,
# )