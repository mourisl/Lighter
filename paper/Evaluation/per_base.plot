set terminal png truecolor enhanced linewidth 2
set terminal jpeg truecolor enhanced linewidth 2
set output "per_base.jpg"

set key bottom 

set xlabel "position"
set ylabel "correct ratio(%)"

set yrange [99:100]

plot "per_base" using 1:2 with line title "Original", "per_base" using 1:3 with line title "Quake", "per_base" using 1:4 with line title "Musket", "per_base" using 1:5 with line title "Bless", "per_base" using 1:6 with line title "Lighter"
