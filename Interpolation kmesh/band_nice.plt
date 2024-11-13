set xtics ( \
    "L"   -0.167407, \
    "A"    0.000000, \
    "L"    0.167407 )
set ytics ( \
    "H"   -0.167407, \
    "A"    0.000000, \
    "H"    0.167407 )

set xrange [-0.167407: 0.167407]
set yrange [-0.167407: 0.167407]

set terminal pdfcairo enhanced font "DejaVu" transparent fontscale 1 size 5.00in, 7.50in
set output "band.pdf"
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set zlabel "E - E_{VBM} (eV)"
set zrange [-0.5: 0.5]
unset key
set ztics 0.1 scale 1 nomirror out
set mytics 2
set terminal qt persist
set parametric

# Set up the 3D view
set view 30, 60
set grid
set style data lines

# Iterate over all partition files to plot on the same 3D graph
splot for [i=1:10] sprintf("band_partition_%d.dat", i) using 1:2:3 with lines lt 1 lw 0.5 title sprintf("Partition %d", i)
