
clear all
close all

f = fopen('../Average_Collision_Angles.dat'); 
d = fscanf(f,'%f %f %f',[3,inf]); 
d = d'; 

figure
h1 = semilogy(d(:,1),d(:,2),'k','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Collision Energy [eV]'); 
ylabel('Scattering Angle Mean [deg]'); 

figure
h2 = semilogy(d(:,1),d(:,3),'g','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Collision Energy [eV]'); 
ylabel('Scattering Angle Standard Deviation [deg]'); 


