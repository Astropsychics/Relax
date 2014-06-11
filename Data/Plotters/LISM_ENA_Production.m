
clear all
close all

f1 = fopen('../lism_production_hp.dat'); 
d1 = fscanf(f1,'%f %f',[2,inf]); 
d1 = d1'; 

r  = d1(:,1); 
f  = d1(:,2); 

AU = 1.5e11; 

figure
h1 = loglog(f,r,'k'); 
set(h1,'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('ENA Production [1/m^3/s]'); 
ylabel('Distance from Star [AU]'); 
print -depsc2 ./Plots/LISM_ENA_Production.eps
print -djpeg100 ./Plots/LISM_ENA_Production.jpeg

dr = r(2) - r(1); 
dr = dr*AU; 

tot = 0; 
for i=1:length(r)
	tot = tot + f(i)*dr; 
end

H_tot = tot; 

for i=1:length(r)
	tot = 0; 
	for j=1:i
		tot = tot + f(j)*dr; 
	end
	Cum_H(i) = tot/H_tot; 
end

figure
h2 = semilogy(Cum_H,r,'g'); 
set(h2,'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Cumulative ENA Production Percentage');
ylabel('Distance from Star [AU]'); 
print -depsc2 ./Plots/LISM_Cumulative_ENA_Production.eps
print -djpeg100 ./Plots/LISM_Cumulative_ENA_Production.jpeg

