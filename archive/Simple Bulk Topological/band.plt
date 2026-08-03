ef=  5.99843566
set xtics ( \
"{/Symbol \107}"   -0.837034, \
"M"   -0.000000, \
"K"    0.483262, \
"{/Symbol \107}"    1.449785, \
"A"    1.907019, \
"L"    2.744053, \
"H"    3.227315, \
"A"    4.193838 )
set xrange [   -0.837034:    4.193838]
set terminal pdfcairo enhanced font "Times New Roman"  transparent fontscale 1 size 4.75in, 4.25in
set output "topologicalband.pdf"
set encoding iso_8859_1
set size square
set ylabel "E-E_{F} (eV)"
set yrange [ -6 : 4.0 ]
#set xrange [ 1.907019 : 4.193838 ]
unset key
set ytics 1.0 scale 1 nomirror out
set mytics 2
set parametric
set trange [-10:10]
plot "band.dat" u 1:($2-ef) with l lt 1 lw 1 lc rgb "blue",\
   -0.000000,t with l lt 2  lc -1,\
    0.483262,t with l lt 2  lc -1,\
    1.449785,t with l lt 2  lc -1,\
    1.907019,t with l lt 2  lc -1,\
    2.744053,t with l lt 2  lc -1,\
    3.227315,t with l lt 2  lc -1,\
t,0 with l lt 2  lc -1
