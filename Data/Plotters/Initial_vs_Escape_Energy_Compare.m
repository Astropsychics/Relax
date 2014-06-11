
clear all
close all

f = fopen('../Initial_Energy.dat');
I = fscanf(f,'%f',[1,inf]);
I = I';

f = fopen('../Escape_Energy_Distribution.dat');
E = fscanf(f,'%f',[1,inf]);
E = E';

NH = 25; 

[yI,xI] = hist(I,NH); 
[yE,xE] = hist(E,NH); 

dxI = xI(2) - xI(1); 
dxE = xE(2) - xE(1); 

CI = 1/(sum(yI)*dxI); 

yI = CI*yI; 
yE = CI*yE; 

figure
plot(xI,yI,'k',xE,yE,'r','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Energy [eV]')
ylabel('Frequency [1/eV]')

fsw = fopen('./Data/Initial_Energy_Dist.dat','w'); 
fes = fopen('./Data/Escape_Energy_Dist.dat','w'); 

for i=1:length(xI)
	fprintf(fsw,'%e %e\n', xI(i), yI(i)); 
	fprintf(fes,'%e %e\n', xE(i), yE(i)); 
end





