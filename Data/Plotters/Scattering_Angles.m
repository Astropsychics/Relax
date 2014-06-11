
clear all
close all

f_co2 = fopen('../planet_3d_CO2_Scatt_Angle.dat'); 
f_o   = fopen('../planet_3d_O_Scatt_Angle.dat'); 
f_he  = fopen('../planet_3d_He_Scatt_Angle.dat'); 
f_h   = fopen('../planet_3d_H_Scatt_Angle.dat'); 

d_co2 = fscanf(f_co2,'%f %f',[2,inf]); 
d_o   = fscanf(f_o,'%f %f',[2,inf]); 
d_he  = fscanf(f_he,'%f %f',[2,inf]); 
d_h   = fscanf(f_h,'%f %f',[2,inf]); 

d_co2 = d_co2'; 
d_o   = d_o'; 
d_he  = d_he'; 
d_h   = d_h'; 

co2_max = max(d_co2(:,1)); 

figure
plot(d_co2(:,1),d_co2(:,2),'k.')
title('He+CO2 Scattering Angles')
xlabel('Lab Frame Collision Energy [eV]')
ylabel('Lab Frame Scattering Angle [deg]')

figure
plot(d_o(:,1),d_o(:,2),'b.')
title('He+O Scattering Angles')
xlabel('Lab Frame Collision Energy [eV]')
ylabel('Lab Frame Scattering Angle [deg]')

figure
plot(d_he(:,1),d_he(:,2),'r.')
title('He+He Scattering Angles')
xlabel('Lab Frame Collision Energy [eV]')
ylabel('Lab Frame Scattering Angle [deg]')

figure
plot(d_h(:,1),d_h(:,2),'g.')
title('He+H Scattering Angles')
xlabel('Lab Frame Collision Energy [eV]')
ylabel('Lab Frame Scattering Angle [deg]')







