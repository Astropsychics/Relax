
clear all
close all

f = fopen('../planet_3d_Atmosphere_Density.dat'); 
d = fscanf(f,'%f %f %f %f %f',[5,inf]); 
d = d'; 

CC = 1e-6; 

z   = d(:,1)/1000; 
co2 = d(:,2)*CC; 
o   = d(:,3)*CC; 
he  = d(:,4)*CC; 
h   = d(:,5)*CC; 

figure
h = semilogx(co2,z,'k',o,z,'b',he,z,'r',h,z,'g'); 
axis([1 10^20 0 1000])
set(h,'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Density [1/cm^3]')
ylabel('Altitude [km]')
legend('CO2', 'O', 'He', 'H','Location','Best')
print -depsc2 ./Plots/Atmosphere_Density.eps
print -djpeg100 ./Plots/Atmosphere_Density.jpeg

