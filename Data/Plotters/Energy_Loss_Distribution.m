
clear all
close all

PRINT_ON = 1; 

f = fopen('../planet_3d_Energy_Loss_Height.dat'); 
d = fscanf(f,'%f %f',[2,inf]); 
d = d'; 

dE = d(:,1); 
z  = d(:,2)/1000; 

zmin = min(z); 
%zmax = max(z); 
zmax = 400; 

Nz   = 51; 
dz   = (zmax-zmin)/(Nz-1); 
Z    = zmin:dz:zmax; 

for i=1:(Nz-1)
	z1 = Z(i); 
	z2 = Z(i+1); 
	c  = 0;
	t  = 0;  
	for j=1:length(z)
		if ( (z1 <= z(j)) && (z(j) <= z2) )
			t = t+dE(j); 
			c = c+1; 
		end
	end
	EE(i)	= t/c; 
	ZZ(i) = z1 + (z2-z1)/2; 
end

figure
h1 = semilogx(dE,z,'k.', EE,ZZ,'r');
set(h1,'LineWidth',2.5) 
set(gca,'FontSize',16)
xlabel('Energy Loss per Collision [eV]','FontSize',16)
ylabel('Alitude [km]','FontSize',16)
if (PRINT_ON == 1)
	print -depsc2 ./Plots/Energy_Lost_Collisions_Averages.eps
	print -djpeg100 ./Plots/Energy_Lost_Collisions_Averages.jpeg
end

figure
h3 = semilogx(EE,ZZ,'r*',EE,ZZ,'k'); 
set(h3,'LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Average Energy Loss per Collision [eV]','FontSize',16)
ylabel('Alitude [km]','FontSize',16)
if (PRINT_ON == 1)
	print -depsc2 ./Plots/Energy_Lost_Averages.eps
	print -djpeg100 ./Plots/Energy_Lost_Averages.jpeg
end


% Write final data to file for better production plots
f1 = fopen('./Grace_Data/Energy_Loss_Height.dat','w'); 
for i=1:length(EE)
	fprintf(f1,'%f %f\n',EE(i),ZZ(i)); 
end

