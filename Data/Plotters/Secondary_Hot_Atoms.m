
clear all
close all

fid = fopen('../planet_3d_Secondary_Hot.dat'); 
dat = fscanf(fid,'%f %f %f %f %f %f',[6,inf]); 
dat = dat'; 

E_proj = dat(:,1); 
E_targ = dat(:,2); 
Height = dat(:,3); 
vx     = dat(:,4); 
vy     = dat(:,5); 
vz     = dat(:,6);

N = length(E_proj); 

figure
[Y1,X1] = hist(E_proj, 100);
Y1 = Y1/N;
h1 = plot(X1,Y1,'b');   
set(h1,'LineWidth',2.5)
set(gca,'FontSize',16)
title('Energy of Projectiles which create Secondary Hot Atoms','FontSize',16)
xlabel('Energy [eV]','FontSize',16)
ylabel('Percent Frequency','FontSize',16)


figure
[Y2,X2] = hist(E_targ, 100); 
Y2 = Y2/N; 
h2 = semilogy(X2,Y2,'b'); 
set(h2,'LineWidth',2.5)
set(gca,'FontSize',16)
title('Energy of Secondary Hot Atoms','FontSize',16)
xlabel('Energy [eV]','FontSize',16)
ylabel('Percent Frequency','FontSize',16)
axis([0 max(E_targ) 1e-6 1])


figure
[Y3,X3] = hist(Height, 100);
Y3 = Y3/N;
X3 = X3/1000;  
h3 = plot(Y3,X3,'b'); 
set(h3,'LineWidth',2.5);
set(gca,'FontSize',16)
title('Height of Creation of Secondary Hot Atoms')
ylabel('Height [km]','FontSize',16)
xlabel('Percent Frequency','FontSize',16)
print -depsc2 ./Plots/Secondary_Hot_Creation_Height.eps

figure
[Y4,X4] = hist(vx, 100); 
Y4 = Y4/N; 
h4 = plot(X4,Y4,'b'); 
set(h4,'LineWidth',2.5); 
set(gca,'FontSize',16)
title('Vertical Velocity Component of Secondary Hot Atoms','FontSize',16)
xlabel('V_{x}','FontSize',16)
ylabel('Percent Frequency','FontSize',16)
axis([-1 1 0 1.01*max(Y4)])

up = 0; 
dn = 0; 
for i=1:length(Y4)
	if (X4(i) > 0)
		up = up + 1; 
	else
		dn = dn + 1; 
	end
end

up = up/length(Y4); 
dn = dn/length(Y4); 
xx = [up dn];
lb1 = sprintf('%2.0f%% Upward', up*100);  
lb2 = sprintf('%2.0f%% Downward', dn*100);  
lb = {lb1,lb2};
ex = [0 1];  
 	
figure
h5 = pie(xx,ex,lb); 
set(h5,'LineWidth',2.5); 
set(h5(2),'fontsize',16,'position',[-0.6 0],'color','k'); 
set(h5(4),'fontsize',16,'position',[0.6 0],'color','K'); 
colormap winter
title('Vertical Flux Component Secondary Hot Atoms','FontSize',16)

 
