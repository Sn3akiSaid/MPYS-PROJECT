# Plot customizations
set xlabel "k (Å^{-1})" font "Times New Roman,10"
set ylabel "E-E_{F} (eV)" font "Times New Roman,10" offset 1.2,0
set border  # Display borders on left and bottom only


# Adjusted Arrows for the Top Axis
set style arrow 1 nohead dashtype 2 lw 1.5  # Dashed arrow style
set arrow from -0.01, graph 1.03 to -0.1, graph 1.03 arrowstyle 1 lc rgb "black" # A -> H
set arrow from  0.01, graph 1.03 to  0.1, graph 1.03 arrowstyle 1 lc rgb "black" # A -> L
 

# Set the style for the filled head (no line, just the head)
set style arrow 2 head filled size screen 0.02,7 lw 1.5# Solid head, no line
# Draw the solid filled heads separately at the ends
set arrow from -0.1, graph 1.03 to -0.12, graph 1.03 nohead arrowstyle 2 lc rgb "black" # Head at H
set arrow from 0.1, graph 1.03 to 0.12, graph 1.03 nohead arrowstyle 2 lc rgb "black" # Head at L

# Ensure the top ticks are displayed correctly
set x2tics ( \
    "L" -0.13, \
    "A"  0.000000, \
    "H"  0.13 ) scale 0.01 offset 0,-0.5 font "Times New Roman" 10

set xrange [-0.14: 0.14]
set x2range [-0.14: 0.14]
set terminal pngcairo enhanced font "Times New Roman"  transparent fontscale 1 size 5.00in, 4.50in dpi 300
set output "interpolated_orbitals10.png"
set encoding iso_8859_1
set size ratio 0 1.0,1.0

set yrange [ -0.5 : 0.5 ]
#set key above Left box 3
set xtics -0.15,0.05,0.15 font "Times New Roman,10"
set ytics -0.4,0.2,0.4 font "Times New Roman,10"
#set mytics 5
set parametric
set trange [-10:10]
plot \
     "band_partition_10.dat" using ($1):($2):($4*0.01) with circles lt 5 lc rgb 'blue' fs transparent solid 0.5 noborder notitle "p_z", \
     "band_partition_10.dat" using ($1):($2):($5*0.01) with circles lt 18 lc rgb 'green' fs transparent solid 0.5 noborder notitle "d_{yz}+d_{xz}", \
     "band_partition_10.dat" using ($1):($2):($3*0.01) with circles lt 2 lc rgb 'red' fs transparent solid 0.5 noborder notitle "p_x+p_y", \
     "band_partition_10.dat" u 1:($2) with l lc "yellow" lt -1 lw 1 notitle ,\
   -0.000000,t with l lt 2  lc -1 notitle,\
t,0 with l lt 2  lc -1 notitle
