# Set terminal for output (optional, can also output to a PNG or PDF file)
set terminal pngcairo enhanced size 800, 600
set output 'contour_plot.png'

# Set labels for the axes
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Z'

# Set the view to 2D (for contour plot)
set view map

# Set the color palette for the contour plot
set palette model RGB defined ( 0 "blue", 1 "green", 2 "yellow", 3 "red")

# Define the grid size and limits (based on your data's range)
set xrange [-0.09:0.1]  # You can adjust these values to match your data range
set yrange [-0.09:0.1]  # Adjust this as well based on your data range

# Interpolate the data into a grid using dgrid3d
set dgrid3d 10,10   # This will create a 20x20 grid (adjust as needed)

# Turn on contour plotting
set contour base  # Enable contour drawing at the base of the graph

# Set the levels of contours (you can adjust this to match your data range)
set cntrparam levels incremental 0.5, 0.5, 3.0  # Start at 0.5, increment by 0.5, max 3.0

# Enable smooth interpolation for the surface to ensure a smooth color transition
set pm3d interpolate 10,10   # This smooths the transition between grid points

# Generate contour plot using the pm3d style
# Generate band contour plot of specified band partition
splot 'band_partition_3.dat' using 1:2:3 with pm3d notitle 
#splot for [i=1:10]sprintf("band_partition_%d.dat", i) using 1:2:3 with pm3d notitle