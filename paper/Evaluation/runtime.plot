set terminal png truecolor enhanced linewidth 3

set output "runtime.png"

set xlabel "# of threads"
set ylabel "time(s)"

plot "quake.runtime" with linespoints title "Quake", "musket.runtime" with linespoints title "Musket", "lighter.runtime" with linespoints title "Lighter"
