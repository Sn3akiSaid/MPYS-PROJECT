ef=  4.18903780
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
set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 7.50in
set output "band.pdf"
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set ylabel "E-E_{VBM} (eV)"
set yrange [ -6 : 4.0 ]
unset key
set ytics 1.0 scale 1 nomirror out
set mytics 2
set parametric
set trange [-10:10]
plot "band.dat" u 1:($2-ef) with l lt 1 lw 3,\
   -0.000000,t with l lt 2  lc -1,\
    0.483262,t with l lt 2  lc -1,\
    1.449785,t with l lt 2  lc -1,\
    1.907019,t with l lt 2  lc -1,\
    2.744053,t with l lt 2  lc -1,\
    3.227315,t with l lt 2  lc -1,\
t,0 with l lt 2  lc -1
