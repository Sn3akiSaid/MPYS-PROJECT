ef=  4.026
set xtics ( \
"K"   -0.966523, \
"{/Symbol \107}"   -0.000000, \
"K"    0.966523 )
set xrange [ -0.167407: 0.167407]
set title "k_z=0.00{/Symbol \160}/c"
set terminal pdfcairo enhanced font "Helvetica"  transparent fontscale 1 size 5.00in, 4.50in
set output "interpolatedweight-kz=0.00.pdf"
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set ylabel "E-E_{F} (eV)"
set yrange [ -4.5 : -3.6 ]
set key above Left box 3
set ytics 0.5 scale 1
set mytics 5
set parametric
set trange [-10:10]
plot \
     "spin_partition_9.dat" using ($1):($2-ef):(abs($3)*0.01) with circles lt 2 lc rgb 'red' fs transparent solid 0.5 noborder notitle "p_x+p_y", \
     "spin_partition_9.dat" using ($1):($2-ef):(abs($4)*0.01) with circles lt 5 lc rgb 'blue' fs transparent solid 0.5 noborder notitle "p_z", \
     "spin_partition_9.dat" using ($1):($2-ef):(abs($5)*0.01) with circles lt 18 lc rgb 'green' fs transparent solid 0.5 noborder notitle "d_{yz}+d_{xz}", \
     "spin_partition_9.dat" u 1:($2-ef) with l lc "yellow" lt -1 lw 1 notitle ,\
   -0.000000,t with l lt 2  lc -1 notitle,\
t,0 with l lt 2  lc -1 notitle
