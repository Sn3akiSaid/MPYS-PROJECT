gap = 'gap.dat'

stats gap using 2 nooutput

print "Min_E:", STATS_min

#Plots the gap energy against alpha values
set title "energy vs alpha"
set xlabel "alpha"
set ylabel "energy"

set terminal qt persist

plot gap using 1:2 with linespoints lt 1 lw 0.5