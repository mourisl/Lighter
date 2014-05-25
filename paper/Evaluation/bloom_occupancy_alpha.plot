set terminal png truecolor enhanced linewidth 2
#set terminal jpeg truecolor enhanced linewidth 2
set terminal postscript eps color enhanced linewidth 2
set output "bloom_occupancy_alpha.eps"

set key top left 

set xlabel "alpha"
set ylabel "occupancy(%)"

set yrange [0:100]

plot "bloom_occupancy_alpha" using 1:2 with linespoint title "35x Table A", \
	"bloom_occupancy_alpha" using 1:3 with linespoint title "35x Table B", \
	"bloom_occupancy_alpha" using 1:4 with linespoint title "70x Table A", \
	"bloom_occupancy_alpha" using 1:5 with linespoint title "70x Table B", \
	"bloom_occupancy_alpha" using 1:6 with linespoint title "140x Table A", \
	"bloom_occupancy_alpha" using 1:7 with linespoint title "140x Table B"
