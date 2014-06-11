
clear all

LY = 1.58e-5; 

f = fopen('../LISM_ENA_Coll_Hist.dat'); 
d = fscanf(f,'%f %f',[2,inf]); 
d = d';

fcx = fopen('../LISM_ENA_CX_Coll_Hist.dat'); 
fat = fopen('../LISM_ENA_Atom_Coll_Hist.dat'); 
dcx = fscanf(fcx,'%f %f',[2,inf]); 
dat = fscanf(fat,'%f %f',[2,inf]); 

dcx = dcx'; 
dat = dat'; 

Hx = d(:,1);
Hy = d(:,2);  

clear d

f = fopen('../LISM_ENA_Coll_Percent.dat'); 
d = fscanf(f,'%f %f',[2,inf]); 
d = d';

Px = d(:,1); 
Py = d(:,2); 

figure
loglog(Hx,Hy,'k','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Number of Collisions Until Thermalization'); 
ylabel('Frequency'); 
print -depsc2 ./Plots/LISM_ENA_Coll_Hist.eps

%figure
%plot(Px,Py,'k','LineWidth',2.5); 
%set(gca,'FontSize',16); 
%xlabel('Number of Collisions Until Thermalization'); 
%ylabel('Percent Frequency'); 
%print -depsc2 ./Plots/LISM_ENA_Coll_Per.eps

NNN = length(dcx(:,1)); 

figure
loglog(dcx(2:NNN,1),dcx(2:NNN,2),'r','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Number of CX Collisions Until Thermalization'); 
ylabel('Frequency'); 
print -depsc2 ./Plots/LISM_ENA_CX_Coll_Hist.eps

figure
loglog(dat(:,1),dat(:,2),'g','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Number of Elastic Atom-Atom Collisions Until Thermalization'); 
ylabel('Frequency'); 
print -depsc2 ./Plots/LISM_ENA_CX_Coll_Hist.eps

figure
loglog(Hx,Hy,'k', dat(:,1), dat(:,2), 'g', dcx(2:NNN,1), dcx(2:NNN,2), 'r', 'LineWidth', 2.5); 
set(gca,'FontSize',16); 
xlabel('Number of Collsions During Transport')
ylabel('Frequency'); 
legend('Total','Atom-Atom','Atom-Ion','Location','Best'); 
print -depsc2 ./Plots/LISM_ENA_Collisions.eps

FILE_WRITE = 1; 
if (FILE_WRITE == 1)
	ftot = fopen('./Plots/Total_Collisions.dat','w'); 
	fat  = fopen('./Plots/Atomic_Collisions.dat','w'); 
	fcx  = fopen('./Plots/CX_Collisions.dat','w'); 
	for i=1:length(Hx)
		fprintf(ftot,'%f %f\n', Hx(i), Hy(i)); 
	end
	for i=1:length(dat(:,1))
		fprintf(fat,'%f %f\n', dat(i,1), dat(i,2)); 
	end
	for i=2:NNN
		fprintf(fcx,'%f %f\n', dcx(i,1), dcx(i,2)); 
	end
end	
