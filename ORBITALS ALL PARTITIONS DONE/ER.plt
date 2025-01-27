i=1
ymin=0
ymax=0.2
file = sprintf("band_partition_%d.dat", i)
set table "filter.dat"
plot file u 1:2:(($2>=ymin&&$2<=ymax) ? $2 : 1/0) w table
unset table
stats "filter.dat" using 2 nooutput
print STATS_min 
