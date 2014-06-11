
clear all

f = fopen('../LISM_ENA_Energy_Dist.dat'); 
d = fscanf(f,'%f %f',[2,inf]); 
d = d'; 

LY = 1.58e-5; 

r = d(:,1)*LY; 
E = d(:,2); 

figure
loglog(r,E,'k','LineWidth',2.5); 
%semilogy(r,E,'k','LineWidth',2.5); 
%plot(r,E,'k','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Distance from Sun [LY]'); 
ylabel('Average Energy [eV]'); 
print -depsc2 ./Plots/LISM_ENA_Energy_Dist.eps


