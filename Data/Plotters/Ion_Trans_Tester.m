
clear all
close all

f1 = fopen('../ion_trans_test_end_r.dat'); 
f2 = fopen('../ion_trans_test_end_v.dat'); 

d1 = fscanf(f1,'%f %f %f',[3,inf]); 
d2 = fscanf(f2,'%f %f %f',[3,inf]); 

d1 = d1'; 
d2 = d2'; 

AU = 1.5e11; % number of [m] in an [AU]
j  = 1; 

for i=1:length(d1)
	R(i,1) = d1(i,1)/AU; 
	R(i,2) = d1(i,2)/AU; 
	R(i,3) = d1(i,3)/AU; 
	r0(i)  = sqrt( R(i,1)^2 + R(i,2)^2 + R(i,3)^2 ); 
	if (r0(i) <= 20)
		Sr0(j) = r0(i); 
		j      = j+1; 
	end
	v0     = sqrt(d2(i,1)^2 + d2(i,2)^2 + d2(i,3)^2); 
	V(i,1) = d2(i,1)/v0; 
	V(i,2) = d2(i,2)/v0; 
	V(i,3) = d2(i,3)/v0; 
end

[syr,sxr] = hist(Sr0,100); 

[yr,xr]   = hist(r0,100); 

[Vxf, Vx] = hist(V(:,1),100); 
[Vyf, Vy] = hist(V(:,2),100); 
[Vzf, Vz] = hist(V(:,3),100); 

N   = length(V(:,1)); 
Vxf = Vxf/N; 
Vyf = Vyf/N; 
Vzf = Vzf/N; 

figure
h1 = plot3(R(:,1),R(:,2),R(:,3),'k.','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('X [AU]')
ylabel('Y [AU]')
zlabel('Z [AU]')

figure
h2 = plot3(V(:,1),V(:,2),V(:,3),'g.','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('vx')
ylabel('vy')
zlabel('vz')

figure
h3 = semilogy(xr,yr,'k*','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Displacement [AU]'); 
ylabel('Frequency'); 

figure
h33 = plot(sxr,syr,'k*','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Displacement [AU]'); 
ylabel('Frequency'); 

figure
h4 = plot(Vx,Vxf,'k',Vy,Vyf,'b',Vz,Vzf,'g','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('Unit Velocity Component'); 
ylabel('Frequency'); 
legend('Vx','Vy','Vz','Location','Best'); 

