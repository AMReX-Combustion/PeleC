set terminal wxt 1
#set terminal png enhanced color font 'Courier,16'
#set output 'N2_NIST_density_60bar.png'

set key top right

set xlabel "Temperature [K]"  
set ylabel "Density [Kg/m3]"

plot "NIST_N2_60bar.txt" u 1:3 title "NIST" w p lc 1 ps 1.5 pt 1,\
     "Test_SRK_vsNIST_N2.txt" u 2:($3*1000.0) title "SRK-PelePhysics" w l lc -1 lt 6 lw 2.0


set terminal wxt 2
#set terminal png enhanced color font 'Courier,16'
#set output 'N2_NIST_Cp_60bar.png'

set key top right

set xlabel "Temperature [K]"  
set ylabel "Specific heat [J/gK]"

plot "NIST_N2_60bar.txt" u 1:9 title "NIST" w p lc 1 ps 1.5 pt 1,\
     "Test_SRK_vsNIST_N2.txt" u 2:($4*1e-7) title "SRK-PelePhysics" w l lc -1 lt 6 lw 2.0


set terminal wxt 3
#set terminal png enhanced color font 'Courier,16'
#set output 'N2_NIST_SpeedofSound_60bar.png'

set key top right

set xlabel "Temperature [K]"  
set ylabel "Speed of Sound [m/s]"

plot "NIST_N2_60bar.txt" u 1:10 title "NIST" w p lc 1 ps 1.5 pt 1,\
     "Test_SRK_vsNIST_N2.txt" u 2:($6*1e-2) title "SRK-PelePhysics" w l lc -1 lt 6 lw 2.0


