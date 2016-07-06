set grid xtics ytics
set yrange [0.6:2.8]
plot "Pressure0.txt" u 1:3  w l lw 1.5 title 'Spring constant = 0',\
     "Pressure100.txt" u 1:3  w l lw 1.5 title 'Spring constant = 100',\
     "Pressure200.txt" u 1:3  w l lw 1.5 title 'Spring constant = 200',\
     "Pressure400.txt" u 1:3  w l lw 1.5 title 'Spring constant = 400',\
     "Pressure800.txt" u 1:3  w l lw 1.5 title 'Spring constant = 800'
set terminal postscript eps color enhanced "Arial" 25
set output "Pressure_group.eps"
replot

set grid xtics ytics
set yrange [0:700]
plot "Pressure0.txt" u 1:2  w l lw 1.5 title 'Spring constant = 0',\
     "Pressure100.txt" u 1:2  w l lw 1.5 title 'Spring constant = 100',\
     "Pressure200.txt" u 1:2  w l lw 1.5 title 'Spring constant = 200',\
     "Pressure400.txt" u 1:2  w l lw 1.5 title 'Spring constant = 400',\
     "Pressure800.txt" u 1:2  w l lw 1.5 title 'Spring constant = 800'
set terminal postscript eps color enhanced "Arial" 25
set output "Volume_group.eps"
replot
