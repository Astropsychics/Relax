
clear all
close all

f0 = fopen('../lism_X_xyz.dat'); 
f1 = fopen('../lism_all_count_XY.dat'); 
f2 = fopen('../lism_all_count_XZ.dat'); 
f3 = fopen('../lism_all_count_YZ.dat'); 

d0 = fscanf(f0,'%f',[1,inf]); 
d1 = fscanf(f1,'%d %d %f', [3,inf]); 
d2 = fscanf(f2,'%d %d %f', [3,inf]); 
d3 = fscanf(f3,'%d %d %f', [3,inf]); 

%AUtoPC = 1; 
AUtoPC = 4.85e-6; 

d0 = d0'*AUtoPC; 
d1 = d1'; 
d2 = d2'; 
d3 = d3'; 


XYZ = d0; 
NH 	= sqrt(length(d1(:,1))); 

small = 0; 

k = 1; 
for i=1:NH
	for j=1:NH
		if (d1(k,3) == 0)
			R_XY(i,j) = small; 	
		else
			R_XY(i,j) = d1(k,3); 
		end
		if (d2(k,3) == 0)
			R_XZ(i,j) = small; 	
		else
			R_XZ(i,j) = d2(k,3); 
		end
		if (d2(k,3) == 0)
			R_YZ(i,j) = small; 	
		else	
			R_YZ(i,j) = d3(k,3); 
		end
		k = k+1; 
	end
end

dx = XYZ(2)-XYZ(1);

C1 = 1/(sum(sum(R_XY(:,:)))*dx*dx);  
C2 = 1/(sum(sum(R_XZ(:,:)))*dx*dx);  
C3 = 1/(sum(sum(R_YZ(:,:)))*dx*dx);  

PL_XY = log10(R_XY*C1); 
PL_XZ = log10(R_XZ*C2); 
PL_YZ = log10(R_YZ*C3); 
P_XY = (R_XY*C1); 
P_XZ = (R_XZ*C2); 
P_YZ = (R_YZ*C3); 

mn   = 0;
rng  = max(max(P_XY(:,:))) - mn;  

NC = 0; 
NN = length(P_XY(:,1)); 
N2 = floor(NN/2)+1; 

v(1:NN,1:NN,1:NN) = 0; 

for i=1:NN
	for j=1:NN
		for k=1:NN
			if (k==N2)
				v(i,j,k) = PL_XY(i,j); 	
			end
			if (i==N2)
				v(i,j,k) = PL_YZ(j,k); 	
			end
			if (j==N2)
				v(i,j,k) = PL_XZ(i,k); 	
			end
		end
	end
end

v(:,:,:)  = 0; 
v(:,:,N2) = PL_XY(:,:); 
v(:,N2,:) = PL_XZ(:,:); 
v(N2,:,:) = PL_YZ(:,:); 

xslice = 0; yslice = 0; zslice = 0; 
figure
slice(XYZ,XYZ,XYZ,v,xslice,yslice,zslice)
%whitebg('k')
set(gca,'FontSize',16)
caxis([-4.5,-0.5])
colormap(jet)
%colormap(flipud(hot))
h = colorbar; 
ylabel(h,'Probability Density [PC^{-2}]','FontSize',16);
%colormap(hot)
view(3);
axis on;
grid on;
%light;
%lighting phong;
%camlight('left');
shading interp;
xlabel('X [PC]'); 
ylabel('Y [PC]'); 
zlabel('Z [PC]'); 
hold on
contour(XYZ,XYZ,PL_XY)


figure
contourf(XYZ,XYZ,P_XY)
colormap(jet)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('X [PC]')
ylabel('Y [PC]')
h = colorbar;
az = 0; el = 60; 
%view(az,el)
ylabel(h,'Probability Density [PC^{-2}]','FontSize',16);
print -depsc2 ./Plots/LISM_Ave_Engy/PD_XY.eps

figure
contourf(XYZ,XYZ,P_XZ)
colormap(jet)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('X [PC]')
ylabel('Z [PC]')
h = colorbar; 
az = 90; el = 60; 
%view(az,el)
ylabel(h,'Probability Density [PC^{-2}]','FontSize',16);
print -depsc2 ./Plots/LISM_Ave_Engy/PD_XZ.eps

figure
contourf(XYZ,XYZ,P_YZ)
colormap(jet)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('Y [PC]')
ylabel('Z [PC]')
h = colorbar; 
az = -90; el = 45; 
%view(az,el)
ylabel(h,'Probability Density [PC^{-2}]','FontSize',16);
print -depsc2 ./Plots/LISM_Ave_Engy/PD_YZ.eps


figure
contourf(XYZ,XYZ,PL_XY)
colormap(hot)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('X [PC]')
ylabel('Y [PC]')
h = colorbar;
az = 0; el = 60; 
%view(az,el)
ylabel(h,'Log10 Probability Density [1/PC^{2}]','FontSize',16);
print -depsc2 ./Plots/LISM_Ave_Engy/PD_Log_XY.eps

figure
contourf(XYZ,XYZ,PL_XZ)
colormap(hot)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('X [PC]')
ylabel('Z [PC]')
h = colorbar; 
az = 90; el = 60; 
%view(az,el)
ylabel(h,'Log10 Probability Density [1/PC^{2}]','FontSize',16);
print -depsc2 ./Plots/LISM_Ave_Engy/PD_Log_XZ.eps

figure
contourf(XYZ,XYZ,PL_YZ)
colormap(hot)
set(gca,'YDir','Normal')
set(gca,'FontSize',16)
xlabel('Y [PC]')
ylabel('Z [PC]')
h = colorbar; 
az = -90; el = 45; 
%view(az,el)
ylabel(h,'Log10 Probability Density [1/PC^{2}]','FontSize',16);
print -depsc2 ./Plots/LISM_Ave_Engy/PD_Log_YZ.eps


