
clear all
close all

f1 = fopen('../ena_production_hp.dat'); 
f2 = fopen('../ena_production_hepp.dat'); 
d1 = fscanf(f1,'%f %f',[2,inf]); 
d2 = fscanf(f2,'%f %f',[2,inf]); 
d1 = d1'; 
d2 = d2'; 

H_z  = d1(:,1); 
H_p  = d1(:,2); 
He_z = d2(:,1); 
He_p = d2(:,2); 

z    = H_z; 
p    = H_p + He_p; 

figure
h = semilogx(H_p,z,'k',He_p,z,'b',p,z,'r--'); 
set(h,'LineWidth',2.5); 
set(gca,'FontSize',16)
axis([1e-2 1e8 100 2000])
xlabel('ENA Production [1/m^3/s]','FontSize',16)
ylabel('Altitude [km]','FontSize',16)
legend('H','He','All ENAs')
title('Mars')
print -depsc2 ./Plots/Mars_ENA_Production.eps
print -djpeg100 ./Plots/Mars_ENA_Production.jpeg

WRITE_FILES = 1; 

if (WRITE_FILES == 1) 
	%% convert to 1/cm^3 for "production" plots
	CCC = 1e-6; 
	fH  = fopen('./Grace_Data/H_P.dat','w'); 
	fHe = fopen('./Grace_Data/He_P.dat','w'); 
	fT  = fopen('./Grace_Data/Tot_P.dat','w'); 
	for i=1:length(z)
		fprintf(fH, '%e %e\n', H_p(i)*CCC, z(i)); 
		fprintf(fHe,'%e %e\n', He_p(i)*CCC, z(i)); 
		fprintf(fT, '%e %e\n', p(i)*CCC, z(i)); 
	end
end

dz = z(2) - z(1);
dz = dz*1000; 

tot = 0; 
for i=1:length(z)
	tot = tot + He_p(i)*dz; 
end

He_C = tot; 

tot = 0; 
for i=1:length(z)
	tot = tot + H_p(i)*dz; 
end

H_C = tot; 

for i=1:length(z)
	tot1 = 0; 
	tot2 = 0; 
	for j=1:i
		tot1 = tot1 + He_p(j)*dz; 
		tot2 = tot2 + H_p(j)*dz; 
	end
	Cum_He(i) = tot1/He_C;
	Cum_H(i)  = tot2/H_C;  
end

figure
h1 = semilogy(Cum_He,z,'b',Cum_H,z,'k'); 
set(h1,'LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Cumulative ENA Production Percentage','FontSize',16)
ylabel('Altitude [km]','FontSize',16)
legend('He','H','Location','NorthWest')
title('Mars')
print -depsc2 ./Plots/Mars_Cumulative_ENA_Production.eps
print -djpeg100 ./Plots/Mars_Cumulative_ENA_Production.jpeg

WRITE_ON = 1; 

if (WRITE_ON == 1)
	fo1 = fopen('../Mars_H_ENA_Start.dat','w'); 
	fo2 = fopen('../Mars_He_ENA_Start.dat','w'); 
	for i=1:length(z)
		fprintf(fo1,'%e %e\n', Cum_H(i), z(i)); 
		fprintf(fo2,'%e %e\n', Cum_He(i), z(i)); 
	end
end


