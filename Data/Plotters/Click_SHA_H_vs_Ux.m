
clear all
close all

fid = fopen('../planet_SHA_height_vs_Ux.dat'); 
dat = fscanf(fid,'%f %f %f',[3,inf]); 
dat = dat'; 

Z   = dat(:,1); 
U   = dat(:,2); 
P   = dat(:,3); 

z0  = U(1); 
i   = 2; 
while( U(i) ~= z0 )
	i = i+1; 
end
N_Hist = i-1; 

Uz  = U(1:N_Hist); 

k = 1; 
for i=1:N_Hist
	for j=1:N_Hist
		PP(i,j) = P(k); 
		k = k+1; 
	end
	ZZ(i) = Z(k-1)/1000; 
end

figure
imagesc(Uz,ZZ,PP)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','Normal')
xlabel('Uz')
ylabel('Altitude [km]')	
colorbar

