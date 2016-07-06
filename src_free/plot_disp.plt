set terminal dumb
plot 'Pressure.txt' u 1:2 title 'Expansion room volume.(cm^3)'
plot 'Pressure.txt' u 1:3 title 'Expansion room pressure.(MPa)'

# Eps format is prefered
set terminal postscript eps color enhanced "Arial" 25
set output "Volume.eps"
set yrange [0:550]
set xtics 60
set ytics 100
set grid xtics ytics
plot "Pressure.txt" u 1:2  w l lw 1.5 title 'Expansion room volume.(cc)'

# Also output a png version for convenience
set terminal png font "/home/bao/fonts_bao/Fonts/arial.ttf" 25 size 800,600
set output "Volume.png"
replot

# Reset the output
set output

set terminal postscript eps color enhanced "Arial" 25
set output "Pressure.eps"
set xtics 60
set ytics 0.2
set yrange [0.6:2.5]
# Reset the grid setting.
unset grid
set grid xtics ytics
plot "Pressure.txt" u 1:3  w l lw 1.5 title 'Expansion room pressure.(MPa)'

set terminal png font "/home/bao/fonts_bao/Fonts/arial.ttf" 25 size 800,600
set output "Pressure.png"
replot

set output
