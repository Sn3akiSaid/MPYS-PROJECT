set terminal wxt persist
set size square
datafile = 'k_points_for_alpha_0.01_perturbed.dat'
# Set axis labels
set xlabel 'k_x (Å^{-1})'
set ylabel 'k_y (Å^{-1})'

# First get the actual data ranges
stats datafile using 1 nooutput
x_min = STATS_min
x_max = STATS_max

stats datafile using 2 nooutput
y_min = STATS_min
y_max = STATS_max

stats datafile using 3 nooutput
z_min = STATS_min
z_max = STATS_max

# Set plot ranges (as specified)
set xrange [-0.1:0.1]
set yrange [-0.1:0.1]

# Set up color palette

set palette defined (z_min "blue", z_max "red")
set cbrange [z_min:z_max]
set cblabel 'k_z (Å^{-1})'

# Enable colorbox and set its position
set colorbox

# Print ranges for debugging
print sprintf("Data ranges: x[%.3f:%.3f] y[%.3f:%.3f] z[%.3f:%.3f]", \
             x_min, x_max, y_min, y_max, z_min, z_max)

# Plot the data
splot datafile using 1:2:3 with points pt 7 ps 1.0 palette title ''