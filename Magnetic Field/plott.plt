#Plot customizations
set xlabel "k_{x} (Å^{-1})" font "Times New Roman,10"
set ylabel "E-E_{F} (eV)" font "Times New Roman,10" offset 1.2,0
set xrange [-0.14: 0.14]
#set x2range [-0.14: 0.14]
set yrange [ -0.5 : 0.5 ]
#stats "band_partition_1.dat" using 3
#set cbrange [-1:1]

set palette defined (-1 "red", 0 "white", 1 "blue")
set style line 1 dashtype 3 pointsize default lw 3
set terminal pdfcairo enhanced font "Times New Roman" fontscale 1 size 5.00in, 4.50in
set output "By,kx,triv-FIXED.pdf"
plot "band_partition_1.dat" u 1:2:3 w l ls 1 lc palette notitle, \
     "band_partition_1.dat" u 1:4:5 w l lw 2 lc palette notitle
