set terminal gif animate delay 5
set output 'Gifs/Ising104_4.gif'
data = 'Data/Data1.dat'
stats data

set xrange [0.5:60.5]
set yrange [0.5:60.5]

do for [bb=0:int(STATS_blocks)-1]{
    set title sprintf('Frame:%03.0f',bb)
    plot data index(bb) pt 20 ps 1.3
}