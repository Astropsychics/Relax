
clear all
close all

fid = fopen('../Rand_Init_Energy.dat'); 
dat = fscanf(fid,'%f %f %f %f',[4,inf]); 
dat = dat'; 

r = dat(:,1); 
y = dat(:,2); 
m = dat(:,3); 
V = dat(:,4); 

[yr,xr] = hist(r,100); 
[yy,xy] = hist(y,100); 
[ym,xm] = hist(m,100); 
[yV,xV] = hist(V,100); 

figure
plot(xr,yr,'k','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('RND'); 
ylabel('Frequency'); 


figure
plot(xy,yy,'b','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('y0'); 
ylabel('Frequency'); 


figure
plot(xm,ym,'c','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('slope'); 
ylabel('Frequency'); 


figure
plot(xV,yV,'g','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Velocity'); 
ylabel('Frequency'); 


