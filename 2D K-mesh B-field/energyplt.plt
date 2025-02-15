#This plot is for the 3D energy 

set terminal wxt persist
#set terminal pdfcairo enhanced font "Times New Roman" transparent fontscale 0.5 size 10.00in, 10.00in
set view 75, 120, 0.5, 1
set output "3D_energy_plot.pdf"
set zrange [ -0.6 : 0.7 ]
set xrange [ -0.4 : 0.4 ]
set yrange [ -0.4 : 0.4 ]
set ztics -0.6, 0.2, 0.6

#set tics scale 0.5
set size square # This makes the plot fill out the page more

set lmargin at screen 0.1
set rmargin at screen 0.9
set bmargin at screen 0.1
set tmargin at screen 0.9


# Define a palette for the z-axis (energy)
#set palette defined ( -0.6 "blue", 0 "white", 0.6 "red" )
# or
set palette rgbformulae 33,13,10
# or
#set palette rgbformulae 21,22,23
# or
#set palette rgbformulae 7,5,15
# or
#set palette cubehelix
# or
#set palette viridis
# Enable pm3d for smooth surfaces
set pm3d interpolate 20,20
set dgrid3d 50,50 qnorm 2
set style data pm3d
#set pm3d depthorder hidden3d # Enable depth ordering and hidden3d for transparency

# Set transparency level for surfaces
set style fill transparent solid 0.3


# Plot the data
datafile = "Energy_part_0.2B1.dat"
#datafile = "Energy_part1.dat"

# Plot the data
splot datafile u 1:2:3 w pm3d notitle,\
      datafile u 1:2:4 w pm3d notitle,\
      datafile u 1:2:5 w pm3d notitle,\
      datafile u 1:2:6 w pm3d notitle

#splot "Energy_part1.dat" u 1:2:3 with points pt 7 ps 1 lc palette notitle,\
#      "Energy_part1.dat" u 1:2:4 with points pt 7 ps 1 lc palette notitle,\
#      "Energy_part1.dat" u 1:2:5 with points pt 7 ps 1 lc palette notitle,\
#      "Energy_part1.dat" u 1:2:6 with points pt 7 ps 1 lc palette notitle

#splot "Energy_part_0.2B1.dat" u 1:2:3 with points pt 7 ps 1 lc palette notitle,\
#      "Energy_part_0.2B1.dat" u 1:2:4 with points pt 7 ps 1 lc palette notitle,\
#      "Energy_part_0.2B1.dat" u 1:2:5 with points pt 7 ps 1 lc palette notitle,\
#      "Energy_part_0.2B1.dat" u 1:2:6 with points pt 7 ps 1 lc palette notitle

