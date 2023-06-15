set xlabel "x"
set ylabel "u"
set grid
set terminal pngcairo size 600, 400
set output "result.png"
plot "result.dat" u 1:2 w l title "Численное решение", "result.dat" u 1:3 w l title "Аналитическое решение"
