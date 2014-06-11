
clear all

N   = 1e5; 
I0  = 1e6; 

fid = fopen('../planet_SHA_Production.dat'); 
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
xlabel('Average Nascent SHA Energy [eV]')
ylabel('Altitude [km]')

figure
plot(dC*I0/N,Z,'k','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Average Nascent SHA Flux [cm^{-2} sec^{-1}]')
ylabel('Altitude [km]')
