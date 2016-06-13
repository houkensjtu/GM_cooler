set terminal dumb
plot 'Pressure.txt' u 1:2 title 'Expansion room volume.(cm^3)'
plot 'Pressure.txt' u 1:3 title 'Expansion room pressure.(MPa)'

set grid xtics ytics
set yrange [0.6:2.5]
plot "Pressure.txt" u 1:3  w l lw 1.5 title 'Expansion room pressure.(MPa)'
set terminal png font "/home/bao/fonts_bao/Fonts/arial.ttf" 25 size 800,600
set output "Pressure.png"
replot

set grid xtics ytics
set autoscale y
plot "Pressure.txt" u 1:2  w l lw 1.5 title 'Expansion room volume.(cc)'
set terminal png font "/home/bao/fonts_bao/Fonts/arial.ttf" 25 size 800,600
set output "Volume.png"
replot
