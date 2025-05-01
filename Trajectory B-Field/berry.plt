# berry_curvature.plt

# Set the terminal type and output file
set terminal wxt persist size 1000,800 enhanced font "Arial,12"
# set size 1000,1000
#set output "berry_curvature.pdf"

# Set the title and labels
set title "Berry Curvature Near Fermi Energy"
set xlabel "k_x"
set ylabel "k_y" offset -1,0
# set zlabel "k_z"
set cblabel "k_z"
# Set the style for vectors
data = "berryinWSM_10.dat"

stats data u 8 nooutput

set xrange [-0.06:0.06]
set yrange [-0.06:0.06]
set zrange[0.441:0.443]
# set zrange[0.44:0.4415]
# Set the viewing angle
set view 0, 360, 1.2, 1.2
set palette defined(-0.001"red",0"white",+0.001"blue")
set cbrange [-0.0005:0.0005]

set style arrow 1 head filled size 0.1,0.2,60 lc palette#lc #rgb "blue"
# Set the key/legend position
# set key outside
# Define energy filtering parameters - replace these with your desired values
energy_minimum = STATS_min  # Set your minimum energy threshold here
dE = 0.05                   # Set your energy window width here
# Main plot command with energy filtering
splot data u 1:2:3:($4/$7*0.003):($5/$7*0.003):($6/$7*0.00):($6/$7*0.001) w vectors arrowstyle 1 not
# splot data u \
#     ($8 > energy_minimum && $8 < energy_minimum+dE ? $1 : 1/0):\
#     ($8 > energy_minimum && $8 < energy_minimum+dE ? $2 : 1/0):\
#     ($8 > energy_minimum && $8 < energy_minimum+dE ? $3 : 1/0):\
#     ($8 > energy_minimum && $8 < energy_minimum+dE ? $4/$7*0.001 : 1/0):\
#     ($8 > energy_minimum && $8 < energy_minimum+dE ? $5/$7*0.001 : 1/0):\
#     ($8 > energy_minimum && $8 < energy_minimum+dE ? $6/$7*0.00 : 1/0):\
#     ($8 > energy_minimum && $8 < energy_minimum+dE ? $7 : 1/0):($6/$7*0.001)\
#     w vectors arrowstyle 1 notitle
# Close the output file
set output
