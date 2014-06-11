
clear all
close all

fid = fopen('../Mars_Density_Test.dat'); 
dat = fscanf(fid,'%f %f %f %f %f %f %f %f %f',[9,inf]); 
dat = dat'; 

C   = 1e-6; 

Z   = dat(:,1); 
H   = dat(:,2)*C; 
He  = dat(:,3)*C; 
O   = dat(:,4)*C; 
Ar  = dat(:,5)*C; 
H2  = dat(:,6)*C; 
N2  = dat(:,7)*C; 
CO  = dat(:,8)*C; 
CO2 = dat(:,9)*C; 

Z   = Z/1000; 

figure
semilogx(H,Z,'k',He,Z,'b',O,Z,'g',Ar,Z,'r', ...
         H2,Z,'k--',N2,Z,'b--',CO,Z,'g--',CO2,Z,'r--', ...
				'LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Density [1/cm^3]')
ylabel('Altitude [km]')
axis([1e0 1e15 min(Z) max(Z)])
legend('H','He','O','Ar','H2','N2','CO','CO2','Location','Best')



