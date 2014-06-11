
clear all
close all

fid = fopen('../planet_3d_Trace.dat'); 
dat = fscanf(fid,'%f %f %f',[3,inf]); 
dat = dat'; 

x = dat(:,1); 
y = dat(:,2); 
z = dat(:,3); 

mkm = 1/1000; 

figure
plot3(x,y,z,'k.', 260, 0, 0, 'go')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
set(gca,'FontSize',16)

figure
plot(y*mkm,x*mkm,'b.', 0, 260, 'go')
xlabel('Y [km]','FontSize',16)
ylabel('X [km]','FontSize',16)
set(gca,'FontSize',16)

figure
plot(z*mkm,x*mkm,'r.', 0, 260, 'go')
xlabel('Z [km]','FontSize',16)
ylabel('X [km]','FontSize',16)
set(gca,'FontSize',16)



