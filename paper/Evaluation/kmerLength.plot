set terminal png truecolor enhanced linewidth 2
set terminal jpeg truecolor enhanced linewidth 2

set output "kmerLength.jpg"

set key bottom 

set xlabel "k-mer length"
set ylabel "percent(%)"

set yrange [0:100]

plot "kmerLength" using 1:2 with linespoint title "Recall", "kmerLength" using 1:3 with linespoint title "Precision", "kmerLength" using 1:4 with linespoint title "F-score", "kmerLength" using 1:5 with linespoint title "Gain"
