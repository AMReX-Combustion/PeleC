set terminal wxt 1
#set terminal png enhanced color font 'Courier,20'
#set output 'N2_dpdT_AnalyticalVsFD.png'

set xlabel "Temperature [K]"  
set ylabel "dpdT |_{tau,Y}"

set format y  "%g"

plot "Test_SRK_1e-6_Derivatives.txt" u 2:5 title "FD - delta = 1e-6" w lp lc rgb "green" lt 2 lw 2.0 ps 0.5 pt 2,\
     "Test_SRK_1e-6_Derivatives.txt" u 2:4 title "Analytical" w l lc rgb "red" lt -1 lw 2.0 



set terminal wxt 2
#set terminal png enhanced color font 'Courier,20'
#set output 'N2_dpdtau_AnalyticalVsFD.png'

set xlabel "Temperature [K]"  
set ylabel "dpdtau |_{T,Y}"

set format y  "%g"

plot "Test_SRK_1e-6_Derivatives.txt" u 2:7 title "FD - delta = 1e-6" w lp lc rgb "green" lt 2 lw 2.0 ps 0.5 pt 2,\
     "Test_SRK_1e-6_Derivatives.txt" u 2:6 title "Analytical" w l lc rgb "red" lt -1 lw 2.0 



set terminal wxt 3
#set terminal png enhanced color font 'Courier,20'
#set output 'N2_dhdtau_AnalyticalVsFD.png'

set xlabel "Temperature [K]"  
set ylabel "dhdtau |_{T,Y}"

set format y  "%g"

plot "Test_SRK_1e-6_Derivatives.txt" u 2:9 title "FD - delta = 1e-6" w lp lc rgb "green" lt 2 lw 2.0 ps 0.5 pt 2,\
     "Test_SRK_1e-6_Derivatives.txt" u 2:8 title "Analytical" w l lc rgb "red" lt -1 lw 2.0 


set terminal wxt 4
#set terminal png enhanced color font 'Courier,20'
#set output 'N2_SumHyk_vs_HM.png'

set xlabel "Temperature [K]"  
set ylabel "Enthalpy [erg/g]"

set format y  "%g"

set key bottom right

plot "Test_SRK_1e-6_Derivatives.txt" u 2:10 title "Analytical - hm" w l lc rgb "red" lt -1 lw 2.0,\
     "Test_SRK_1e-6_Derivatives.txt" u 2:11 title "Sum (Yk hk)" w lp lc rgb "green" lt 2 lw 2.0 ps 0.5 pt 2


set terminal wxt 5
#set terminal png enhanced color font 'Courier,20'
#set output 'N2_damdT_AnalyticalVsFD.png'

set xlabel "Temperature [K]"  
set ylabel "damdT "

set format y  "%g"

set key bottom right

plot "Test_SRK_1e-6_Derivatives.txt" u 2:13 title "FD - delta = 1e-6" w lp lc rgb "green" lt 2 lw 2.0 ps 0.5 pt 2,\
     "Test_SRK_1e-6_Derivatives.txt" u 2:12 title "Analytical" w l lc rgb "red" lt -1 lw 2.0


set terminal wxt 6
#set terminal png enhanced color font 'Courier,20'
#set output 'N2_d2amdT2_AnalyticalVsFD.png'

set xlabel "Temperature [K]"  
set ylabel "d^{2} am / dT^{2} "

set format y  "%g"

set key bottom right

plot "Test_SRK_1e-6_Derivatives.txt" u 2:15 title "FD - delta = 1e-6" w lp lc rgb "green" lt 2 lw 2.0 ps 0.5 pt 2,\
     "Test_SRK_1e-6_Derivatives.txt" u 2:14 title "Analytical" w l lc rgb "red" lt -1 lw 2.0
     				    
