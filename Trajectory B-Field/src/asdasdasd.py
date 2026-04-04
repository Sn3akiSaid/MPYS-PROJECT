import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev
from mpl_toolkits.mplot3d import Axes3D

# Load data
data = np.loadtxt("trysm.dat")
x, y, z = data[:, 0], data[:, 1], data[:, 2]

# Handle stacking by taking the mean z for duplicate (x, y) pairs
from collections import defaultdict

xyz_dict = defaultdict(list)
for i in range(len(x)):
    xyz_dict[(x[i], y[i])].append(z[i])

# Compute the representative z (mean, median, or first occurrence)
x_unique, y_unique, z_unique = [], [], []
for (x_val, y_val), z_vals in xyz_dict.items():
    x_unique.append(x_val)
    y_unique.append(y_val)
    z_unique.append(np.mean(z_vals))  # You can use np.median(z_vals) instead

x, y, z = np.array(x_unique), np.array(y_unique), np.array(z_unique)

# Sort points (optional: based on angle in polar coordinates)
angles = np.arctan2(y - np.mean(y), x - np.mean(x))
sorted_indices = np.argsort(angles)
x, y, z = x[sorted_indices], y[sorted_indices], z[sorted_indices]

# Ensure closure by repeating the first point at the end
x = np.append(x, x[0])
y = np.append(y, y[0])
z = np.append(z, z[0])

# Fit a periodic B-spline
tck, u = splprep([x, y, z], s=0, per=True)

# Evaluate the spline with more points
u_fine = np.linspace(0, 1, 200)
x_smooth, y_smooth, z_smooth = splev(u_fine, tck)

# Plot the result
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, 'ro', label="Processed Data Points")
ax.plot(x_smooth, y_smooth, z_smooth, 'b-', label="Fitted Loop")
ax.legend()
plt.show()
