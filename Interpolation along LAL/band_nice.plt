set xtics ( \
"H"   -0.167407, \
"A"    0.000000, \
"L"    0.167407 )
set xrange [   -0.167407:    0.167407]
set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 7.50in
set output "band.pdf"
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set ylabel "E-E_{VBM} (eV)"
set yrange [ -0.5 : 0.5 ]
unset key
set ytics 0.1 scale 1 nomirror out
set mytics 2
set terminal qt persist
set parametric
set trange [-10:10]

#iterates over all the partition files to plot on the same graph
plot for [i=1:10] sprintf("band_partition_%d.dat", i) u 1:2 with l lt 1 lw 0.5, \
    0.0, t with l lt 2 lc -1, \
    t, 0.0 with l lt 2 lc -1
