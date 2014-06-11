
clear all

fEL 	= fopen('../LISM_ENA_dE_NColl_EL_H.dat'); 
fCX 	= fopen('../LISM_ENA_dE_NColl_CX_H.dat'); 
dELH 	= fscanf(fEL,'%f %f',[2,inf]); 
dCXH 	= fscanf(fCX,'%f %f',[2,inf]); 
dELH 	= dELH'; 
dCXH 	= dCXH'; 
fEL 	= fopen('../LISM_ENA_dE_NColl_EL_He.dat'); 
fCX 	= fopen('../LISM_ENA_dE_NColl_CX_He.dat'); 
dELHe = fscanf(fEL,'%f %f',[2,inf]); 
dCXHe = fscanf(fCX,'%f %f',[2,inf]); 
dELHe = dELHe'; 
dCXHe = dCXHe'; 

BOTH = 0; 

if (BOTH == 1)
figure
loglog(dEL(:,1),dEL(:,2),'k','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Number of Collisions')
ylabel('Average Energy [eV]')
title('Elastic Collisions')

figure
loglog(dCX(:,1),dCX(:,2),'k','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Number of Collisions')
ylabel('Average Energy [eV]')
title('Charge Exchange Collisions')
end

figure
loglog(dELH(:,1),dELH(:,2),'b',dCXH(:,1),dCXH(:,2),'g', ...
dELH(:,1),(dELH(:,2)+dCXH(:,2)),'r--','LineWidth',2.5)
set(gca,'FontSize',16)
title('Hydrogen Projectiles')
xlabel('Number of Collisions')
ylabel('Average Energy [eV]')
legend('Elastic Collisions','Charge-Exchange Collisions','Total','Location','Best')
print -depsc2 ./Plots/Average_Energy_NColl_H.eps

figure
loglog(dELHe(:,1),dELHe(:,2),'b',dCXHe(:,1),dCXHe(:,2),'g', ...
dELHe(:,1),(dELHe(:,2)+dCXHe(:,2)),'r--','LineWidth',2.5)
set(gca,'FontSize',16)
title('Helium Projectiles')
xlabel('Number of Collisions')
ylabel('Average Energy [eV]')
legend('Elastic Collisions','Charge-Exchange Collisions','Total','Location','Best')
print -depsc2 ./Plots/Average_Energy_NColl_He.eps

tot   = dELH(:,2)+dELHe(:,2)+dCXH(:,2)+dCXHe(:,2); 
totH  = dELH(:,2)+dCXH(:,2); 
totHe = dELHe(:,2)+dCXHe(:,2); 

figure
loglog(dELHe(:,1),totH,'b',dCXHe(:,1),totHe,'g', ...
dELHe(:,1),tot,'r--','LineWidth',2.5)
set(gca,'FontSize',16)
title('All Projectiles')
xlabel('Number of Collisions')
ylabel('Average Energy [eV]')
legend('All Hydrogen','All Helium','Total','Location','Best')
print -depsc2 ./Plots/Average_Energy_NColl.eps



