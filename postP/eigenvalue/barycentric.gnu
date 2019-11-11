# using gnuplot 5.2
set style line 1 lt 1 lw 3.0 ps 1.25 pt 9 palette
set style line 10 lt 1 lw 3.0 ps 1.25 pt 9 linecolor rgb "black"


unset label
set terminal pdf enhanced
set output 'barycentric.pdf'
set size 0.8,0.8
set size ratio 0.866
set yrange [0.0:0.866]
set xrange [0.0:1.0]
set cbrange [-0.005:0.005]
set noborder
set noxtics
set noytics
unset xlabel
unset ylabel
set cblabel 'r [m]'
set arrow 1 from 0.0,0.0 to 0.5,(sqrt(3.0)/2.0) nohead front lt -1 lw 2.0
set arrow 2 from 0.0,0.0 to 1.0,0.0             nohead front lt -1 lw 2.0
set arrow 3 from 1.0,0.0 to 0.5,(sqrt(3.0)/2.0) nohead front lt -1 lw 2.0
set label 1 '{/:Bold x}_{2c}' at 0.0-0.035,0.0-0.02
set label 2 '{/:Bold x}_{1c}' at 1.0+0.01,0.0-0.02
set label 3 '{/:Bold x}_{3c}' at 0.5-0.02,(sqrt(3.0)/2.0)+0.026

# using bX bY and color choice, e.g., temperature, coords, etc
plot 'barycentric.dat' using ($2):($3):($1) with points ls 1 ps .5 notitle
