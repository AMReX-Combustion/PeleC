set key font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set key spacing 2.0
set xrange [0:1]
plot 'pelec_soln.dat' u 1:($2/9.6e-4) w l lw 3 title "PeleC",\
'Lv_Ihme_JCP_2014' u 1:2 w p ps 2 pt 5 title "Lv and Ihme, JCP,2014 (N2 left)",\
'Grogan_Ihme_ShockWaves_2020' u 1:2 w p ps 2 pt 7 title "Grogan and Ihme, Shock Waves, 2020 (HE left)"
