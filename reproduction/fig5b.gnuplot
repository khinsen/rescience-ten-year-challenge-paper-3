set terminal pngcairo
set output "fig5b.png"
set key bmargin
set yrange [0.:1.]
plot "anisotropy_acoustic_1IEE_2_2_2.plot" w l dt 5 lc rgb "black", "anisotropy_acoustic_1IEE_3_3_3.plot" w l dt 4 lc rgb "black", "anisotropy_acoustic_1IEE_4_4_4.plot" w l dt 3 lc rgb "black", "anisotropy_acoustic_1IEE_5_5_5.plot" w l dt 2 lc rgb "black", "anisotropy_acoustic_1IEE_10_10_10.plot" w l dt 5 lc rgb "black", "anisotropy_acoustic_1IEE_20_20_20.plot" w l dt 1 lc rgb "black"
