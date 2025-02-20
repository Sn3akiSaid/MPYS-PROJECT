#This plot is for the 3D energy 

set terminal wxt persist
#set terminal pdfcairo enhanced font "Times New Roman" transparent fontscale 0.5 size 14.00in, 10.00in
set view 75, 120, 0.5, 1

#set xrange [ -0.5 : 0.5 ]
#set yrange [ -0.5 : 0.5 ]
#set zrange [ -2.5 : 2.5 ]
unset colorbox
# Define the coordinate frame (arrows)
#    - Adjust the "to" coordinates to make them longer or shorter.
# set arrow 1 from 0.2,0.5,0 to 0.3,0.5,0 head filled size screen 0.01,10 lw 2 lc "blue"
# set arrow 2 from 0.2,0.5,0 to 0.2,0.6,0 head filled size screen 0.01,10 lw 2 lc "blue"
# set arrow 3 from 0.2,0.5,0 to 0.2,0.5,0.1 head filled size screen 0.01,10 lw 2 lc "blue"

# Add labels for axes
#set label 1 "k_{x}" at 0.34, 0.5, 0 center
#set label 2 "k_{y}" at 0.2, 0.64, 0 center
#set label 3 "k_{z}"at 0.2, 0.5, 0.12 center
#set tics scale 0.5
#set size square # This makes the plot fill out the page more
# set lmargin at screen -0.2
# set rmargin at screen 0.8
# set bmargin at screen 0.1
# set tmargin at screen 0.95


# Define a palette for the z-axis (energy)
#set palette defined ( -0.6 "blue", 0 "white", 0.6 "red" )
# or
#set palette rgbformulae 33,13,10 #This is the best one
# or
#set palette rgbformulae 21,22,23
# or
#set palette rgbformulae 7,5,15
# or
#set palette cubehelix
# or
#set palette viridis
# Enable pm3d for smooth surfaces
set pm3d #interpolate 2,2
#set dgrid3d 50,50 qnorm 7
set style data pm3d
#set pm3d depthorder hidden3d # Enable depth ordering and hidden3d for transparency

# Set transparency level for surfaces
set style fill transparent solid 0.32

#Output data
set output "Fermi_Surface_Unperturbed_Trial.pdf"
#set output "3D_energy_plot_perturbed_test3.pdf"

# Plot the data
#datafile = "Energy_part_0.2B1.dat"
#datafile = "Energy_part1.dat"
datafile = "perturbed_fermi_energies.dat"
threshold = 1.0
# Plot the data
splot datafile u (($4 >= 0.9 && $4 <= 1) ? $1 : 1/0):(($4 >= 0.9 && $4 <= 1) ? $2 : 1/0):(($4 >= 0.9 && $4 <= 1) ? $3 : 1/0) with points pt 7 ps 0.1 notitle
      #datafile u 1:2:3 w pm3d notitle,\
      #datafile u 1:2:4 w pm3d notitle,\
      #datafile u 1:2:5 w pm3d notitle,\
      #datafile u 1:2:6 w pm3d notitle,\
     

#splot "Energy_part_0.2B1.dat" u 1:2:3 with points pt 7 ps 1 lc palette notitle,\
#      "Energy_part_0.2B1.dat" u 1:2:4 with points pt 7 ps 1 lc palette notitle,\
#      "Energy_part_0.2B1.dat" u 1:2:5 with points pt 7 ps 1 lc palette notitle,\
#      "Energy_part_0.2B1.dat" u 1:2:6 with points pt 7 ps 1 lc palette notitle

