set terminal qt enhanced font "DejaVu" persist
set title "3D Mesh Plot of Band Structure in k-Space"
set xlabel "k_x"
set ylabel "k_y"
set zlabel "Energy (eV)"
set grid
set pm3d
set view map
splot "band_partition_1.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_2.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_3.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_4.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_5.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_6.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_7.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_8.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_9.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_10.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_11.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_12.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_13.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_14.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_15.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_16.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_17.dat" u 1:2:3 with lines palette notitle, \
splot "band_partition_18.dat" u 1:2:3 with lines palette notitle, \
pause -1
