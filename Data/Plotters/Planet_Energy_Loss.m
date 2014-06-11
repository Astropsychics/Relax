
clear all

fid = fopen('../planet_ENA_Energy_Loss.dat'); 
dat = fscanf(fid,'%e %e %e',[3,inf]); 
dat = dat'; 

Z   = dat(:,1)/1000; 
dE  = dat(:,2); 
dC  = dat(:,3); 

for i=1:length(dE)
	if ( dC(i) == 0) 
		ddE(i) = 0; 
	else
		ddE(i) = dE(i)/dC(i); 
	end
end

figure
plot(ddE,Z,'k','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Average ENA Energy Loss per Collision [eV]')
ylabel('Altitude [km]')

PRINT_FILE = 1;
if (PRINT_FILE == 1)
  fo = fopen('../ENA_dE_vs_Z.dat','w');
  for i=1:length(Z)
    fprintf(fo,'%e %e\n',ddE(i), Z(i));
  end
end

