set terminal pngcairo
set output "fig5a.png"
set key bmargin
set yrange [0.:0.0025]
plot "bacoustic_1IEE_2_2_2.plot" using (100*$1) w l dt 5 lc rgb "black", "bacoustic_1IEE_3_3_3.plot" using (100*$1) w l dt 4 lc rgb "black", "bacoustic_1IEE_4_4_4.plot" using (100*$1) w l dt 3 lc rgb "black", "bacoustic_1IEE_5_5_5.plot" using (100*$1) w l dt 2 lc rgb "black", "bacoustic_1IEE_10_10_10.plot" using (100*$1) w l dt 5 lc rgb "black", "bacoustic_1IEE_20_20_20.plot" using (100*$1) w l dt 1 lc rgb "black"
