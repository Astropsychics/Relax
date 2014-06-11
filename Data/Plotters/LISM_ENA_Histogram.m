
clear all
close all

fid = fopen('../LISM_ENA_Relax.dat'); 
dat = fscanf(fid,'%f %f %f %f',[4,inf]); 
dat = dat'; 

x = dat(:,1); 
y = dat(:,2); 
z = dat(:,3); 
C = dat(:,4); 

for i=1:length(dat(:,1))
	r(i) = sqrt(x(i)^2 + y(i)^2 + z(i)^2); 
end

[yr, xr] = hist(r,100); 
[yC, xC] = hist(C,100); 

figure
plot(xr, yr, 'k.')


figure
plot(xC, yC, 'r.')

