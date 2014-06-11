
clear all
close all

fid = fopen('../planet_phase_space.dat'); 
dat = fscanf(fid,'%d %f %f',[3,inf]); 
dat = dat'; 

N   = dat(:,1); 
E   = dat(:,2); 			% [eV]
Z   = dat(:,3)/1000; 	% [km]

% find number of particles

part = 0; 

for i=1:length(N)
	if ( N(i) == 1 ) 
		part = part + 1; 
	end
end

fprintf('%d particles shown in Energy-Altitude phase space\n', part)

figure
plot(E,Z,'k.','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Energy [eV]')
ylabel('Altitude [km]')



