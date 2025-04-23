# berry_curvature.plt

# Set the terminal type and output file
set terminal pngcairo size 1000,800 enhanced font "Arial,12"
set output "berry_curvature.png"

# Set the title and labels
set title "Berry Curvature Near Fermi Energy"
set xlabel "k_x"
set ylabel "k_y"
set zlabel "k_z"

# Set the style for vectors
# set style arrow 1 head filled size 0.08,15,45 lc rgb "blue"
set xrange [-0.05:0.05]
set yrange [-0.05:0.05]

# Set the viewing angle
set view 0, 360, 1.2, 1.2

# Set the key/legend position
set key outside

# Main plot command with energy filtering
splot "curvature.dat" u 1:2:(abs($3-0.432235)<0.0001 ? $3 : 1/0):($6/$9*0.005):($7/$9*0.005):(0) w vectors arrowstyle 1 title "Berry Curvature"

# Close the output file
set output