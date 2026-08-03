

# Set titles and labels
set title "Energy Surface E(x) and E(y)"
set xlabel "kx"
set ylabel "ky"
set zlabel "Energy (E)"

# Set grid
set grid

# Set color palette for surface plot

# Set view angle
set view 180, 0

# Adjust axis ranges based on stats

# Plot the data as a 3D surface
splot for [i=1:10] sprintf("band_partition_%d.dat", i) using 1:2:3 with lines lt 1 lw 0.5 title sprintf("Partition %d", i)
