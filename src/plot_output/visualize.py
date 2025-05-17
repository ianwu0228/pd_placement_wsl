import numpy as np
import matplotlib.pyplot as plt

# Load gradient + position data
data = np.loadtxt("grad_vectors.txt")
x, y = data[:, 0], data[:, 1]
dx, dy = data[:, 2], data[:, 3]
mag = data[:, 4]

# Optional: normalize gradient directions for visualization
norm = np.sqrt(dx**2 + dy**2)
dx_unit = dx / (norm + 1e-8)
dy_unit = dy / (norm + 1e-8)

# === Create a 2D density map by binning cells ===
# Define bin size based on coordinate range
num_bins = 200  # change to 400 or 800 if needed
x_min, x_max = np.min(x), np.max(x)
y_min, y_max = np.min(y), np.max(y)

# Bin for density heatmap
heatmap, xedges, yedges = np.histogram2d(x, y, bins=num_bins, range=[[x_min, x_max], [y_min, y_max]])

# Plot
plt.figure(figsize=(12, 10))

# Plot density heatmap
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.imshow(heatmap.T, origin='lower', cmap='hot', extent=extent, aspect='auto')
plt.colorbar(label='Cell Density (from bin count)')

# Plot cell positions
plt.scatter(x, y, s=4, color='cyan', label='Cells')

# Plot gradient vectors
quiv = plt.quiver(x, y, dx_unit, dy_unit, mag, cmap='cool', scale=30, width=0.002, label='Gradient')
plt.colorbar(quiv, label='Gradient Magnitude')

# Formatting
plt.title("Gradient Field with Cell Density and Positions")
plt.xlabel("X")
plt.ylabel("Y")
plt.axis('equal')
plt.legend()
plt.tight_layout()
plt.savefig("full_visualization.png", dpi=300)
# plt.show()
