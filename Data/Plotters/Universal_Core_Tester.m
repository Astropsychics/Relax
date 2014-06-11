
clear all
close all

fid = fopen('../Uni_test.dat'); 
dat = fscanf(fid,'%f %f %f %f %f %f %f %f %f',[9,inf]); 
dat = dat'; 

E    = dat(:,1); 
U_T  = dat(:,2); 
Q_T  = dat(:,3); 
U_dE = dat(:,4); 
Q_dE = dat(:,5); 
U_TC = dat(:,6); 
Q_TC = dat(:,7); 
rat  = dat(:,8); 
prat = dat(:,9); 

figure
h1 = semilogy(E,U_T,'k.',E,Q_T,'r.'); 
set(h1,'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Collision Energy [eV]')
ylabel('Average Scattering Angle [deg]')
legend('Universal','Quantum')

figure
h2 = semilogy(E,U_dE,'k.',E,Q_dE,'r.'); 
set(h2,'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Collision Energy [eV]')
ylabel('Average Enegy Lost per Collision [eV]')
legend('Universal','Quantum')

figure
h3 = semilogy(E,U_TC,'k.',E,Q_TC,'r.'); 
set(h3,'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Collision Energy [eV]')
ylabel('Total Cross Section [a0^2]')
legend('Universal','Quantum')

figure
h4 = plot(E,rat,'g.'); 
set(h4,'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Collision Energy [eV]')
ylabel('Ratio Q/U TCS*dE')


