#This plot is for the 3D energy 

set terminal wxt persist
#set terminal pdfcairo enhanced font "Times New Roman" transparent fontscale 0.5 size 14.00in, 10.00in
set view 75, 120, 0.5, 1
set size square
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


#set pm3d depthorder hidden3d # Enable depth ordering and hidden3d for transparency

# Set transparency level for surfaces
set style fill transparent solid 0.32

#Output data
set output "Fermi_Surface_Unperturbed_Trial.pdf"

datafile = "perturbed_fermi_energies.dat"
threshold = 1.0
# Plot the data
#splot datafile u (($4 >= 0.5 && $4 <= 0.9) ? $1 : 1/0):(($4 >= 0.5 && $4 <= 0.9) ? $2 : 1/0):(($4 >= 0.5 && $4 <= 0.9) ? $3 : 1/0) with points pt 7 ps 0.1 notitle
splot datafile u 1:2:3:4 w isosurface
     