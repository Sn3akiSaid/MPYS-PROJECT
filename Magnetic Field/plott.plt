#Plot customizations
set xlabel "k (Å^{-1})" font "Times New Roman,10"
set ylabel "E-E_{F} (eV)" font "Times New Roman,10" offset 1.2,0
set xrange [-0.14: 0.14]
#set x2range [-0.14: 0.14]
set yrange [ -0.5 : 0.5 ]

set terminal pdfcairo enhanced font "Times New Roman" fontscale 1 size 5.00in, 4.50in
set output "band1.pdf"
plot "band_partition_1.dat" u 1:2 w l lc rgb "red" lw 1 notitle, \
     "band_partition_1.dat" u 1:3 w l lc rgb "blue" lw 1 notitle
