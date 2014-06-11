
clear all

fid = fopen('../planet_SHA_Energy_Dist.dat'); 
dat = fscanf(fid,'%e %e %e',[3,inf]); 
dat = dat'; 

Z = dat(:,1)/1000; 
E = dat(:,2); 
C = dat(:,3); 

for i=1:length(E)
	if (C(i) == 0)
		dE(i) = 0; 
	else
		dE(i) = E(i)/C(i); 
	end
end

figure
plot(dE,Z,'k','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Average SHA Energy [eV]')
ylabel('Altitude [km]')

FILE_PRINT = 1; 
if (FILE_PRINT == 1)
  fout = fopen('../SHA_Avg_E_vs_Z.dat','w'); 
  for i=1:length(Z)
    fprintf(fout,'%e %e\n',dE(i), Z(i)); 
  end
end


