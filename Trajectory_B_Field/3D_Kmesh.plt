#set terminal pdfcairo enhanced font "Times New Roman"  transparent fontscale 1 size 5.00in, 5.00in
set terminal wxt persist

# Define the base file name and the partition number
#base_file = "mesh_partition_"
#partition_number = 2

# Construct the file name dynamically
#file = sprintf("%s%d.dat", base_file, partition_number)
file = "Spinmesh_TI_3001.dat" 
stats file using 3 nooutput


#set palette defined (0 "black", 0.001 "cyan", 0.015 "green" ,0.02 "yellow", 0.03 "orange", 0.15 "red", 0.2 "white")  # Customize as needed
#set pm3d interpolate 0,0
#set contour base
#set view map
#unset surface
#set cntrparam levels discrete 0,0.1,0.2
#set dgrid3d 100,100

#set table 'contours.dat'
#splot file u 1:2:3 w lines notitle
#unset table
#unset contour
unset cblabel

# Define the color palette for the heat map (energy isocontours)
set size square
# Analyze the 13th column to find its min and max values
# stats outputs "STATS_min" and "STATS_max" for the selected column

# Adjusted Arrows for the Top Axis (x2)
set style arrow 1 nohead dashtype 2 lw 1.5  # Dashed arrow style for x2
set arrow from -0.01, graph 1.04 to -0.08, graph 1.04 arrowstyle 1 lc rgb "black" # A -> H
set arrow from  0.01, graph 1.04 to  0.08, graph 1.04 arrowstyle 1 lc rgb "black" # A -> L

# Set the style for the filled head (no line, just the head for x2)
set style arrow 2 head filled size screen 0.02,7 lw 1.5  # Solid head, no line for x2
# Draw the solid filled heads at the ends for x2
set arrow from -0.07, graph 1.04 to -0.08, graph 1.04 nohead arrowstyle 2 lc rgb "black" # Head at H
set arrow from 0.07, graph 1.04 to 0.08, graph 1.04 nohead arrowstyle 2 lc rgb "black" # Head at L

# Adjusted Arrows for the Right Axis (y2)
# Dashed arrow style for y2
set arrow from graph 1.045, 0.55 to graph 1.045, 0.87 arrowstyle 1 lc rgb "black" # A -> H
set arrow from graph 1.045, 0.1 to graph 1.045, 0.45 arrowstyle 1 lc rgb "black" # A -> L

# Set the style for the filled head (no line, just the head for y2)
# Draw the solid filled heads at the ends for y2
set arrow from graph 1.045, 0.87 to graph 1.045, 0.90 nohead arrowstyle 2 lc rgb "black" # Head at H
set arrow from graph 1.045, 0.13 to graph 1.045, 0.1 nohead arrowstyle 2 lc rgb "black" # Head at L



# Ensure x2 and y2 axes are displayed and linked to the primary axes
#set xrange[0.03:0.05]
#set yrange[0.03:0.05]

set xrange[-0.01:0.01]
set yrange[-0.01:0.01]
# Use the stats results to set the color range dynamically
#set cbrange [0:0.2]
set zrange[0.18:0.2]

#Colorbox positioning
set encoding iso_8859_1
set parametric
set colorbox horizontal user origin 0.23,0.91 # Position it closer to the top (x, y coordinates in the plot space)
set colorbox size 0.57, 0.02 
set cbtics (\
    0, \
    "200" 0.2 ) font "Times New Roman,10" offset 0.22,2.2

set cblabel "E-E_{CBM} (meV)" offset 0.22,4 font "Times New Roman,10"

# Configure axis labels
set xlabel "k_{x} (Å^{-1})" font "Times New Roman,10"
set ylabel "k_{y} (Å^{-1})" font "Times New Roman,10"
set xtics -0.1,0.05,0.1 font "Times New Roman,10"
set ytics -0.1,0.05,0.1 font "Times New Roman,10"

# Set tics for x2 axis
set x2tics ( \
    "L̅" -0.09, \
    "A"  0.000000, \
    "L"  0.09 ) scale 0.01 offset 0,-0.5 font "Times New Roman, 10"

# Set tics for y2 axis
set y2tics ( \
    "H̅" -0.09, \
    "A"  0.000000, \
    "H" 0.09 ) scale 0.01 offset -0.5,0 font "Times New Roman, 10"

# Plot the heatmap and overlay the spin vectors
#'contours.dat' using 1:2:3 lc palette w l not 
#file every 15:15 using 1:2:($14/($17*30)):($15/($17*30)):($16/($17*30)) with vectors head filled size screen 0.02,15,60 lt rgb "black" notitle
#file every 15:15 using 1:2:($4/($7*30)):($5/($7*30)):($6/($7*30)) with vectors head filled size screen 0.02,15,60 lt rgb "black" notitle
#set output "mesh.pdf"
set view 
# Now, plot the heatmap and contours in two separate plot commands
# First, plot the image:
splot for [i=1:20] sprintf("Spinmesh_TI_300%d.dat", i) using 1:2:3 w lines notitle