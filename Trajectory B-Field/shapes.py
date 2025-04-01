import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import curve_fit

# Create a figure and 3D axis
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Load your actual data
data = np.loadtxt("k_points_for_alpha_0.01.dat")
kx = data[:, 0]
ky = data[:, 1]
kz = data[:, 2]
alpha = data[:, 3]

# Calculate actual data ranges
x_min, x_max = np.min(kx), np.max(kx)
y_min, y_max = np.min(ky), np.max(ky)
z_min, z_max = np.min(kz), np.max(kz)
alpha_min, alpha_max = np.min(alpha), np.max(alpha)

print(f"X range: {x_min:.6f} to {x_max:.6f}, span: {x_max-x_min:.6f}")
print(f"Y range: {y_min:.6f} to {y_max:.6f}, span: {y_max-y_min:.6f}")
print(f"Z range: {z_min:.6f} to {z_max:.6f}, span: {z_max-z_min:.6f}")
print(f"Alpha range: {alpha_min:.6f} to {alpha_max:.6f}")

# Find phase by fitting the sine wave to z vs. angle data
x_center = (x_max + x_min) / 2
y_center = (y_max + y_min) / 2

# Calculate the angular position of each data point in the xy plane
angles = np.arctan2(ky - y_center, kx - x_center)
# Ensure angles are in [0, 2π] range
angles = np.mod(angles, 2*np.pi)

# Sort data by angle for better fitting and visualization
sort_idx = np.argsort(angles)
angles_sorted = angles[sort_idx]
z_sorted = kz[sort_idx]

# Function to fit: z = z_center + amplitude * sin(3*t + phase_offset)
def sine_model(t, z_center, amplitude, phase_offset):
    return z_center + amplitude * np.sin(3*t + phase_offset)

# Fit the function to the data
try:
    params, _ = curve_fit(sine_model, angles_sorted, z_sorted)
    z_center_fit, z_amplitude_fit, phase_offset_fit = params
    print(f"Fitted parameters: z_center = {z_center_fit:.6f}, amplitude = {z_amplitude_fit:.6f}, phase = {phase_offset_fit:.6f} rad")
except Exception as e:
    print(f"Fitting failed: {e}")
    # Use estimated parameters if fitting fails
    z_center_fit = (z_max + z_min) / 2
    z_amplitude_fit = (z_max - z_min) / 2
    phase_offset_fit = 0
    print(f"Using estimated parameters: z_center = {z_center_fit:.6f}, amplitude = {z_amplitude_fit:.6f}, phase = {phase_offset_fit:.6f} rad")

# Generate helix model curve using data ranges and fitted phase
radius_x = (x_max - x_min) / 2
radius_y = (y_max - y_min) / 2

t = np.linspace(0, 2*np.pi, 1000)  # Parameter for helix
x_model = x_center + radius_x * np.cos(t)  # X coordinates centered on data
y_model = y_center + radius_y * np.sin(t)  # Y coordinates centered on data
z_model = z_center_fit + z_amplitude_fit * np.sin(3*t + phase_offset_fit)  # Z with fitted phase

# Plot the 3D model curve
ax.plot(x_model, y_model, z_model, 'r-', linewidth=2, alpha=0.7, label='Phase-Adjusted Model')

# Plot the data points
scatter = ax.scatter(kx, ky, kz, 
                    c=alpha,  # Use alpha column for color
                    cmap=cm.plasma,  # Color map
                    s=30,  # Size of points
                    alpha=0.8,  # Opacity of points
                    label='Data')

# Add a color bar
cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
cbar.set_label('Alpha Value')

# Set the axis labels
ax.set_xlabel('kx')
ax.set_ylabel('ky')
ax.set_zlabel('kz')

# Set axis limits with a small margin (5%)
margin = 0.05
x_range = x_max - x_min
y_range = y_max - y_min
z_range = z_max - z_min

ax.set_xlim(x_min - margin * x_range, x_max + margin * x_range)
ax.set_ylim(y_min - margin * y_range, y_max + margin * y_range)
ax.set_zlim(z_min - margin * z_range, z_max + margin * z_range)

# Add legend and title
ax.legend()
plt.title('3D Helix with Phase-Adjusted Model')

# Create 2D projections to visualize the phase adjustment
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot angle vs z to show the phase adjustment
ax1.scatter(angles_sorted, z_sorted, c='blue', s=15, alpha=0.7, label='Data')
t_fine = np.linspace(0, 2*np.pi, 200)
ax1.plot(t_fine, sine_model(t_fine, z_center_fit, z_amplitude_fit, phase_offset_fit), 
         'r-', linewidth=2, label='Fitted Curve')
ax1.set_xlabel('Angle (radians)')
ax1.set_ylabel('z')
ax1.set_title('z vs Angle (Phase Adjustment)')
ax1.legend()
ax1.set_xlim(0, 2*np.pi)

# kx-kz projection
scatter2 = ax2.scatter(kx, kz, c=alpha, cmap=cm.plasma, s=20, alpha=0.8)
# Project model onto kx-kz plane
ax2.plot(x_model, z_model, 'r-', linewidth=2, alpha=0.7)
ax2.set_xlabel('kx')
ax2.set_ylabel('kz')
ax2.set_xlim(x_min - margin * x_range, x_max + margin * x_range)
ax2.set_ylim(z_min - margin * z_range, z_max + margin * z_range)
ax2.set_title('kx-kz Projection with Phase-Adjusted Model')
fig2.colorbar(scatter2, ax=ax2)

plt.tight_layout()
plt.show()