set encoding iso_8859_1
set term pngcairo enhanced color font "Arial-Bold,28" size 1000,1000;

### Plot time in microseconds
set title "System-Name Virtual Dihedral Angles"
set xlabel "Time ({\181}s)"
set ylabel "VDA ({\186})"
set key outside right Left reverse width 2 height 1 font "Arial-Bold,24" maxrows 6
#set xrange [0:5]
#set yrange [-180:180]
#set xtics (0, 1, 2, 3, 4, 5) border nomirror out;
#set ytics (-180, -90, 0, 90, 180) border nomirror out;

set linetype 1 lc rgb "light-coral" 
set linetype 2 lc rgb "dark-cyan" 
set linetype 3 lc rgb "dark-violet"

set output "vda_system-name.png";
plot "system-name_vda_rep1.dat" u ($1/10):($2) w lines t "Rep 1" lw 6, \
"system-name_vda_rep2.dat" u ($1/10):($2) w lines t "Rep 2" lw 6, \
"system-name_vda_rep3.dat" u ($1/10):($2) w lines t "Rep 3" lw 6, \
