set encoding iso_8859_1
set term pngcairo enhanced color font "Arial-Bold,28" size 1000,1000;

set title "Ala1 - Ala2 Distance"
set key inside right Left reverse width 2 height 1 font "Arial-Bold,24" maxrows 6
set xlabel "Time (ns)"
set ylabel "C{/Symbol a} Distance ({\305})"
#set xrange [0:100]
#set yrange [0:10]
#set xtics (0, 20, 40, 60, 80, 100) border nomirror out;
#set ytics (0, 2, 4, 6, 8, 10) border nomirror out;

set linetype 1 lc rgb "light-coral" 
set linetype 2 lc rgb "dark-cyan" 

set output "A1-A2_dist.png";
plot "system_WT_A1-A2_dist.dat" w lines t "System (WT)" lw 6, \
"system_mut_A1-A2_dist.dat" w lines t "System (Mut)" lw 6;
