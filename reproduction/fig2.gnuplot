set terminal png
set output "fig2.png"
set multiplot
set size 1.,0.5
set origin 0.,0.5
set yrange [0.:0.001]
plot "btotal_3LZT_1_1_1.plot" w l
set size 1.,0.5
set origin 0.,0.
set yrange [0.:1.]
plot "anisotropy_total_3LZT_1_1_1.plot" w l
unset multiplot
