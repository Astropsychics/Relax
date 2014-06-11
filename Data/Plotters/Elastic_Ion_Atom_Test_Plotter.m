
clear all; 
close all; 


f1 = fopen('../He_Hp_test.dat'); 
f2 = fopen('../He_Hep_test.dat'); 

d1 = fscanf(f1,'%f %f %f',[3,inf]); 
d2 = fscanf(f2,'%f %f %f',[3,inf]); 

d1 = d1'; 
d2 = d2';

E1 = d1(:,1); 
E2 = d2(:,1); 
T1 = d1(:,2); 
T2 = d2(:,2); 
A1 = d1(:,3); 
A2 = d2(:,3); 

cm = 1e4; 

figure
loglog(E1,T1*cm,'k',E2,T2*cm,'g','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Energy [eV]')
ylabel('Total Cross Section [cm^2]')
legend('He+H^+','He+He^+','Location','Best')

figure
semilogx(E1,A1,'k--',E2,A2,'g--','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Energy [eV]')
ylabel('Average Scattering Angle [deg]')
legend('He+H^+','He+He^+','Location','Best')




