# berry_curvature.plt

# Set the terminal type and output file
set terminal qt persist # size 1000,800 enhanced font "Arial,12"
# set output "berry_curvature.png"

# Set the title and labels
set title "Berry Curvature Near Fermi Energy"
set xlabel "k_x"
set ylabel "k_y"
# set zlabel "k_z"

# Set the style for vectors
# set style arrow 1 head filled size 0.08,15,45 lc rgb "blue"
set xrange [*:*]
set yrange [*:*]
set zrange[0.44:0.4415]
# Set the viewing angle
set view 0, 360, 1.2, 1.2

# Set the key/legend position
# set key outside
data = "berryinWSM_10.dat"
stats data u 8 nooutput
# Define energy filtering parameters - replace these with your desired values
energy_minimum = STATS_min  # Set your minimum energy threshold here
dE = 0.05              # Set your energy window width here
# Main plot command with energy filtering
# splot data u 1:2:3:($4/$7*0.001):($5/$7*0.001):($6/$7*0.00):7 w vector not
splot data u \
    ($8 > energy_minimum && $8 < energy_minimum+dE ? $1 : 1/0):\
    ($8 > energy_minimum && $8 < energy_minimum+dE ? $2 : 1/0):\
    ($8 > energy_minimum && $8 < energy_minimum+dE ? $3 : 1/0):\
    ($8 > energy_minimum && $8 < energy_minimum+dE ? $4/$7*0.001 : 1/0):\
    ($8 > energy_minimum && $8 < energy_minimum+dE ? $5/$7*0.001 : 1/0):\
    ($8 > energy_minimum && $8 < energy_minimum+dE ? $6/$7*0.00 : 1/0):\
    ($8 > energy_minimum && $8 < energy_minimum+dE ? $7 : 1/0) \
    w vector notitle
# Close the output file
set output