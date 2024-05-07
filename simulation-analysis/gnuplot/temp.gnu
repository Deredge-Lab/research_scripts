set encoding iso_8859_1
set term pngcairo enhanced color font "Arial-Bold,28" size 1000,1000;

set title "System-Name C{/Symbol a} Temp (ST)"
set xlabel "# Steps"
set ylabel "Temperature (K)"
set key inside left Left reverse width 2 height 1 font "Arial-Bold,24" maxrows 6
#set xrange [0:50000000]
#set yrange [300:500]
set format x "%.1t{\327}10^%T"
set xtics add ("0" 0) border nomirror out;
set ytics border nomirror out;

set output "system-name_ST_300-500K.png";
plot "system-name_ST_rep1.dat" w lines t "Rep 1" lw 6, \
"system-name_ST_rep2.dat" w lines t "Rep 2" lw 6, \
"system-name_ST_rep3.dat" w lines t "Rep 3" lw 6, \
