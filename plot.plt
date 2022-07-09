set terminal gif animate delay 5
set output 'Gifs/Ising40.gif'
data = 'DataRelI2/Data40.dat'
stats data nooutput

set xrange [0.5:60.5]
set yrange [0.5:60.5]

do for [b=0:19]{
    plot data index(b) pt 20 ps 1.3
}