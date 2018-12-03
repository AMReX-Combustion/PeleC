clear all
close all
format long

%fextract.Linux.gfortran.exe -p plt_2d_002200 -s toto.txt -d 2 -v 'temp y_velocity density X(H2) X(H) X(O) X(O2) X(OH) X(H2O) X(HO2) X(H2O2) X(C) X(CH) X(CH2) X(CH2(S)) X(CH3) X(CH4) X(CO) X(CO2) X(HCO) X(CH2O) X(CH2OH) X(CH3O) X(CH3OH) X(C2H) X(C2H2) X(C2H3) X(C2H4) X(C2H5) X(C2H6) X(HCCO) X(CH2CO) X(HCCOH) X(N) X(NH) X(NH2) X(NH3) X(NNH) X(NO) X(NO2) X(N2O) X(HNO) X(CN) X(HCN) X(H2CN) X(HCNN) X(HCNO) X(HOCN) X(HNCO) X(NCO) X(N2) X(AR) X(C3H7) X(C3H8) X(CH2CHO) X(CH3CHO)'

file_name = "PMF_FROM_PELELM_METHANE_1atm_300K_Phi_0_7.dat"

clear read_data
read_data = dlmread (file_name)(3:end,1:end);

read_data(:,1) = read_data(:,1) *100;
read_data(:,3) = read_data(:,3) *100;
read_data(:,4) = read_data(:,4) *0.001;

textfile = fopen(file_name);
line1 = fgetl(textfile);
line2 = fgetl(textfile);
fclose(textfile);

output_file = "PMF_converted.dat";
fid = fopen (output_file, "w");
fputs (fid, line1);
fputs (fid, "\n");
fputs (fid, line2);
fputs (fid, "\n");
fclose (fid);

dlmwrite(output_file,read_data,"-append","delimiter", " ");