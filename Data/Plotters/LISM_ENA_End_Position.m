
clear all
close all

f = fopen('../LISM_ENA_End_Pos.dat'); 
d = fscanf(f,'%f %f %f',[3,inf]); 
d = d'; 

LY = 1.58e-5; 

x = d(:,1)*LY; 
y = d(:,2)*LY; 
z = d(:,3)*LY; 

N = length(d(:,1)); 

clear d

f = fopen('../LISM_ENA_Diss_Hist.dat'); 
d = fscanf(f,'%f %f',[2,inf]); 
d = d';

Hx = d(:,1)*LY; 
Hy = d(:,2);  

clear d

f = fopen('../LISM_ENA_Diss_Percent.dat'); 
d = fscanf(f,'%f %f',[2,inf]); 
d = d';

Px = d(:,1)*LY; 
Py = d(:,2); 

tit = sprintf('Relaxation Distance\n%d Monte Carlo Particles', N); 

figure
plot3(x,y,z,'g.',0,0,0,'r*','LineWidth',2.5); 
set(gca,'FontSize',16); 
grid on
xlabel('X [LY]'); 
ylabel('Y [LY]'); 
zlabel('Z [LY]');
title(tit);  
print -depsc2 ./Plots/LISM_ENA_End_Position.eps

figure
plot(Hx,Hy,'k','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('X [LY]'); 
ylabel('Frequency'); 
title(tit);  
print -depsc2 ./Plots/LISM_ENA_End_Hist.eps

figure
plot(Px,Py,'k','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('X [LY]'); 
ylabel('Percent Frequency'); 
title(tit);  
print -depsc2 ./Plots/LISM_ENA_End_Per.eps


