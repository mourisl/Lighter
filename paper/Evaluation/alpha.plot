set terminal png truecolor enhanced linewidth 2
#set terminal jpeg truecolor enhanced linewidth 2

set output "alpha.png"

set key bottom 

set xlabel "alpha"
set ylabel "percent(%)"

set yrange [1:100]

plot "alpha" using 1:2 with linespoint title "Recall", "alpha" using 1:3 with linespoint title "Precision", "alpha" using 1:4 with linespoint title "F-score", "alpha" using 1:5 with linespoint title "Gain"

set output "alpha0.png"
set yrange [85:100]

set key bottom center

plot "alpha" using 1:2 with linespoint title "Recall", "alpha" using 1:3 with linespoint title "Precision", "alpha" using 1:4 with linespoint title "F-score", "alpha" using 1:5 with linespoint title "Gain"
