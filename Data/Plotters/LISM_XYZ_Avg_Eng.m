
clear all
close all

f0 = fopen('../lism_X_xyz.dat'); 
f1 = fopen('../lism_all_avg_energy_XY.dat'); 
f2 = fopen('../lism_all_avg_energy_XZ.dat'); 
f3 = fopen('../lism_all_avg_energy_YZ.dat'); 

d0 = fscanf(f0,'%f',[1,inf]); 
d1 = fscanf(f1,'%d %d %f', [3,inf]); 
d2 = fscanf(f2,'%d %d %f', [3,inf]); 
d3 = fscanf(f3,'%d %d %f', [3,inf]); 

AUtoPC = 4.85e-6; 

d0 = d0'*AUtoPC; 
d1 = d1'; 
d2 = d2'; 
d3 = d3'; 


XYZ = d0; 
NH 	= sqrt(length(d1(:,1))); 

k = 1; 
for i=1:NH
	for j=1:NH
		P_XY(i,j) = d1(k,3); 
		P_XZ(i,j) = d2(k,3); 
		P_YZ(i,j) = d3(k,3); 
		k = k+1; 
	end
end

figure
contourf(XYZ,XYZ,P_XY)
colormap(hot)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('X [PC]')
ylabel('Y [PC]')
h = colorbar;
%caxis([0,150])
ylabel(h,'Average Energy [eV]','FontSize',16);
%caxis([-3,2])
print -depsc2 ./Plots/LISM_Ave_Engy/Engy_XY.eps

figure
contourf(XYZ,XYZ,P_XZ)
colormap(hot)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('X [PC]')
ylabel('Z [PC]')
h = colorbar; 
%caxis([0,150])
ylabel(h,'Average Energy [eV]','FontSize',16);
print -depsc2 ./Plots/LISM_Ave_Engy/Engy_XZ.eps

figure
contourf(XYZ,XYZ,P_YZ)
colormap(hot)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('Y [PC]')
ylabel('Z [PC]')
h = colorbar; 
%caxis([0,150])
ylabel(h,'Average Energy [eV]','FontSize',16);
print -depsc2 ./Plots/LISM_Ave_Engy/Engy_YZ.eps


figure
contourf(XYZ,XYZ,log10(P_XY))
colormap(hot)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('X [PC]')
ylabel('Y [PC]')
h = colorbar;
ylabel(h,'Log10 Average Energy [eV]','FontSize',16);
%caxis([-3,2])
print -depsc2 ./Plots/LISM_Ave_Engy/Engy_LOG_XY.eps

figure
contourf(XYZ,XYZ,log10(P_XZ))
colormap(hot)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('X [PC]')
ylabel('Z [PC]')
h = colorbar; 
ylabel(h,'Log10 Average Energy [eV]','FontSize',16);
print -depsc2 ./Plots/LISM_Ave_Engy/Engy_LOG_XZ.eps

figure
contourf(XYZ,XYZ,log10(P_YZ))
colormap(hot)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('Y [PC]')
ylabel('Z [PC]')
h = colorbar; 
ylabel(h,'Log10 Average Energy [eV]','FontSize',16);
print -depsc2 ./Plots/LISM_Ave_Engy/Engy_LOG_YZ.eps

