set terminal png
set output "fig3.png"
set multiplot
set size 1.,0.5
set origin 0.,0.5
set yrange [0.:0.1]
plot "btotal_1IEE_1_1_1.plot" using (100*$1) w l
set size 1.,0.5
set origin 0.,0.
set yrange [0.:1.]
plot "anisotropy_total_1IEE_1_1_1.plot" w l
unset multiplot
