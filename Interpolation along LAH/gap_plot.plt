set terminal pdfcairo enhanced font "Times New Roman"  transparent fontscale 1 size 5.00in, 5.00in

set output "Energy Gap.pdf"
set view
set size square
gap = 'gap.dat'
set style line 1 linecolor rgb "blue" linetype 1 linewidth 1.75 pointtype 7 pointsize 0.75
set style line 2 linecolor rgb "#808080" linetype 2 linewidth 1.75 
stats gap using 2 nooutput

#print "Min_E:", STATS_min
set xrange [0.65517241:0.89655173]
#Plots the gap energy against alpha values
#set title "energy vs alpha"
set xlabel "{/Symbol a}"
set ylabel "E_{Gap} (eV)"

set ytics 0,0.01,0.05 nomirror
set xtics 0.6,0.04,1 nomirror

topweyl=0.79310346 
weyltriv=0.75862068

set label 1 "TI" at 0.84, graph 0.8 center font "Times New Roman, 8"
set label 2 "WSM" at (topweyl+(weyltriv-topweyl)/2), graph 0.8 center font "Times New Roman, 8"
set label 3 "Trivial\nInsulator" at (weyltriv+(0.65517241-weyltriv)/2), graph 0.8 center font "Times New Roman, 8"

set arrow from topweyl, graph 0 to topweyl, graph 1 nohead linestyle 2
set arrow from weyltriv, graph 0 to weyltriv, graph 1 nohead linestyle 2

plot gap using 1:2 w linespoints linestyle 1 notitle
