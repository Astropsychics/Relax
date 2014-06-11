
clear all

JtoeV = 6.2e18; 
M     = 1.7e-27; 

fid = fopen('../LISM_ENA_Start_Energy_Distance.dat'); 
dat = fscanf(fid,'%f %f',[2,inf]); 
dat = dat'; 

ft  = fopen('../../Tables/Raw_SW_Velocity.dat'); 
dt  = fscanf(ft,'%f %f',[2,inf]); 
dt  = dt'; 

dt(:,1) = dt(:,1)/100; 
for i=1:length(dt(:,1))
	Et(i) = 0.5*M*(dt(i,2)*1000)^2; 
end
Et = Et*JtoeV; 

N = length(dat(:,1)); 

[Ey,Ex] = hist(dat(:,1), 100); 
[ry,rx] = hist(dat(:,2), 100); 

Ey = 100*Ey/N; 
ry = 100*ry/N; 

for i=1:length(Ex)
	tot = 0; 
	for j=1:i
		tot = tot + Ey(j); 
	end
	CEy(i) = tot; 
end

for i=1:length(rx)
	tot = 0; 
	for j=1:i
		tot = tot + ry(j); 
	end
	Cry(i) = tot; 
end

figure
h1 = semilogy(Ex,Ey,'r*',Ex,Ey,'k','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Initial ENA Energy [eV]')
ylabel('Percent of Ensemble')
print -depsc2 ./Plots/LISM_Init_ENA_Energy.eps

figure
h2 = semilogy(rx,ry,'r*',rx,ry,'k','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Initial ENA Distance [AU]')
ylabel('Percent of Ensemble')
print -depsc2 ./Plots/LISM_Init_Position.eps

figure
h3 = plot(Ex,CEy,'g', 'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Initial Cumulative ENA Energy [eV]')
ylabel('Percent of Ensemble')
print -depsc2 ./Plots/LISM_Init_Cum_ENA_Energy.eps

figure
h4 = plot(rx,Cry,'b','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Initial Cumulative ENA Location [AU]')
ylabel('Percent of Ensemble')
print -depsc2 ./Plots/LISM_Init_Cum_ENA_Position.eps


