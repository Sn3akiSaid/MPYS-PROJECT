set xtics ( \
"L"   -0.837034, \
"A"    0.000000, \
"L"    0.837034 )
set xrange [   -0.1:    0.1]
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "blue" behind
set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 7.50in
set output "band.pdf"
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set ylabel "E-E_{VBM} (eV)"
set yrange [ -1.5 : -0.5 ]
unset key
set ytics 1.0 scale 1 nomirror out
set mytics 2
set terminal qt persist
set parametric
set trange [-10:10]
plot for [i=2:20] sprintf("band_partition_%d.dat", i) u 1:2 with l lt 1 lw 0.5 lc rgb "white",\

    0.000000,t with l lt 2  lc -1,\
t,0 with l lt 2  lc -1
