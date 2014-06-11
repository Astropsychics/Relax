
clear all
close all

f1 = fopen('../ena_production_hp.dat'); 
f2 = fopen('../ena_production_hepp.dat'); 
d1 = fscanf(f1,'%f %f',[2,inf]); 
d1 = d1'; 
d2 = fscanf(f2,'%f %f',[2,inf]); 
d2 = d2'; 

Z1 = d1(:,1); 
F1 = d1(:,2); 
Z2 = d2(:,1); 
F2 = d2(:,2); 
N = length(Z1); 


figure
plot(F1,Z1,'k',F2,Z2,'g','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Frequency')
ylabel('Altitude [km]')
legend('H','He','Location','Best')



