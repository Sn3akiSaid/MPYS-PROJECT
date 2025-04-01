set terminal wxt persist
set size square
datafile = 'trysm.dat'
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

stats datafile using 4 nooutput
annihilation_min = 0.77
# z_min = annihilation_min + 0.001
annihilation_max = 0.81
# z_max = annihilation_max - 0.001
creation = (STATS_max+STATS_min)/2
# Set plot ranges (as specified)
set xrange [-0.05:0.05]
set yrange [-0.05:0.05]
set zrange [*:*]

# Set up color palette

set palette defined (annihilation_min "blue", creation "white", annihilation_max "red")
set cbrange [0.77:0.807]
set cblabel 'k_z (Å^{-1})'

# Enable colorbox and set its position
set colorbox

# Print ranges for debugging
print sprintf("Data ranges: x[%.3f:%.3f] y[%.3f:%.3f] z[%.3f:%.3f]", \
             x_min, x_max, y_min, y_max, annihilation_min,annihilation_max)

# Plot the data
splot datafile using 1:2:3:4 with points pt 7 ps 1.0 palette title ''