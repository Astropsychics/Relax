
clear all
close all

f1 = fopen('../planet_SHA_height_vs_energy_H.dat'); 
d1 = fscanf(f1,'%f %f %f',[3,inf]); 
d1 = d1'; 
f2 = fopen('../planet_SHA_height_vs_energy_He.dat'); 
d2 = fscanf(f2,'%f %f %f',[3,inf]); 
d2 = d2'; 
f3 = fopen('../planet_SHA_height_vs_energy_O.dat'); 
d3 = fscanf(f3,'%f %f %f',[3,inf]); 
d3 = d3'; 
f4 = fopen('../planet_SHA_height_vs_energy_Ar.dat'); 
d4 = fscanf(f4,'%f %f %f',[3,inf]); 
d4 = d4'; 
f5 = fopen('../planet_SHA_height_vs_energy_H2.dat'); 
d5 = fscanf(f5,'%f %f %f',[3,inf]); 
d5 = d5'; 
f6 = fopen('../planet_SHA_height_vs_energy_N2.dat'); 
d6 = fscanf(f6,'%f %f %f',[3,inf]); 
d6 = d6'; 
f7 = fopen('../planet_SHA_height_vs_energy_CO.dat'); 
d7 = fscanf(f7,'%f %f %f',[3,inf]); 
d7 = d7'; 
f8 = fopen('../planet_SHA_height_vs_energy_CO2.dat'); 
d8 = fscanf(f8,'%f %f %f',[3,inf]); 
d8 = d8'; 

Z   = d1(:,1); 
E   = d1(:,2); 
P1  = d1(:,3); 
P2  = d2(:,3); 
P3  = d3(:,3); 
P4  = d4(:,3); 
P5  = d5(:,3); 
P6  = d6(:,3); 
P7  = d7(:,3); 
P8  = d8(:,3); 

E0  = E(1); 
i   = 2; 
while( E(i) ~= E0 )
	i = i+1; 
end
N_Hist = i-1; 

EE  = E(1:N_Hist); 

k = 1; 
for i=1:N_Hist
	for j=1:N_Hist
		P_H(i,j)   = P1(k); 
		P_He(i,j)  = P2(k); 
		P_O(i,j)   = P3(k); 
		P_Ar(i,j)  = P4(k); 
		P_H2(i,j)  = P5(k); 
		P_N2(i,j)  = P6(k); 
		P_CO(i,j)  = P7(k); 
		P_CO2(i,j) = P8(k); 
		k = k+1; 
	end
	ZZ(i) = Z(k-1)/1000; 
end

for i=1:N_Hist
	C_H(i)   = sum(P_H(:,i)); 
	C_He(i)  = sum(P_He(:,i)); 
	C_O(i)   = sum(P_O(:,i)); 
	C_Ar(i)  = sum(P_Ar(:,i)); 
	C_H2(i)  = sum(P_H2(:,i)); 
	C_N2(i)  = sum(P_N2(:,i)); 
	C_CO(i)  = sum(P_CO(:,i)); 
	C_CO2(i) = sum(P_CO2(:,i)); 
	TOT(i)   = C_H(i) + C_He(i) + C_O(i) + C_Ar(i) + ...
             C_H2(i) + C_N2(i) + C_CO(i) + C_CO2(i); 
end

figure
semilogx(C_H,ZZ,'k',C_He,ZZ,'b', ...
     C_O,ZZ,'g',C_Ar,ZZ,'r', ...
     C_H2,ZZ,'k--',C_N2,ZZ,'b--', ...
     C_CO,ZZ,'g--',C_CO2,ZZ,'r--', ...
		 TOT,ZZ,'c', ...
		'LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('SHA Production')
ylabel('Altitude [km]')
legend('H','He','O','Ar','H2','N2','CO','CO2','Total','Location','Best')
print -djpeg100 ./Plots/SHA_Height_Production.jpeg
 

figure
imagesc(EE,ZZ,P_H)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','Normal')
xlabel('Energy [eV]')
ylabel('Altitude [km]')	
title('SHA H')
colorbar
print -djpeg100 ./Plots/SHA_Height_vs_Energy_H.jpeg

figure
imagesc(EE,ZZ,P_He)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','Normal')
xlabel('Energy [eV]')
ylabel('Altitude [km]')	
title('SHA He')
colorbar
print -djpeg100 ./Plots/SHA_Height_vs_Energy_He.jpeg

figure
imagesc(EE,ZZ,P_O)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','Normal')
xlabel('Energy [eV]')
ylabel('Altitude [km]')	
title('SHA O')
colorbar
print -djpeg100 ./Plots/SHA_Height_vs_Energy_O.jpeg

figure
imagesc(EE,ZZ,P_Ar)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','Normal')
xlabel('Energy [eV]')
ylabel('Altitude [km]')	
title('SHA Ar')
colorbar
print -djpeg100 ./Plots/SHA_Height_vs_Energy_Ar.jpeg

figure
imagesc(EE,ZZ,P_H2)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','Normal')
xlabel('Energy [eV]')
ylabel('Altitude [km]')	
title('SHA H2')
colorbar
print -djpeg100 ./Plots/SHA_Height_vs_Energy_H2.jpeg

figure
imagesc(EE,ZZ,P_N2)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','Normal')
xlabel('Energy [eV]')
ylabel('Altitude [km]')	
title('SHA N2')
colorbar
print -djpeg100 ./Plots/SHA_Height_vs_Energy_N2.jpeg

figure
imagesc(EE,ZZ,P_CO)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','Normal')
xlabel('Energy [eV]')
ylabel('Altitude [km]')	
title('SHA CO')
colorbar
print -djpeg100 ./Plots/SHA_Height_vs_Energy_CO.jpeg

figure
imagesc(EE,ZZ,P_CO2)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','Normal')
xlabel('Energy [eV]')
ylabel('Altitude [km]')	
title('SHA CO2')
colorbar
print -djpeg100 ./Plots/SHA_Height_vs_Energy_CO2.jpeg


