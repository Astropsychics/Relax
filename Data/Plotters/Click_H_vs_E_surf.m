
clear all
close all

fx  = fopen('../planet_X_energy.dat'); 
fy  = fopen('../planet_X_height.dat'); 
RX  = fscanf(fx,'%f',[1,inf]); 
RY  = fscanf(fy,'%f',[1,inf]); 
RX  = RX'; 
RY  = RY'; 
RY  = RY/1000; 

fid = fopen('../planet_height_vs_energy.dat'); 
dat = fscanf(fid,'%d %d %d %f',[4,inf]); 
dat = dat'; 

ck  = dat(:,1); 
Ni  = dat(:,2); 
Nj  = dat(:,3); 
P   = dat(:,4); 

N   = length(P); 

NI  = max(Ni); 
NJ  = max(Nj); 
NC  = N/(NI*NJ); 
CS  = max(ck)/NC;

k = 1; 
figure('units','normalized','outerposition',[0 0 1 1])
for click=1:NC
	tot = 0; 
	for i=1:NI
		for j=1:NJ
			PROB(i,j) = 100*P(k); 
			tot = tot + P(k); 
			k = k+1; 
		end
	end


	tit = sprintf('%04d CLICK\t%5.2f Percent Ensemble in Frame',click*CS,tot*100); 
	suptitle(tit); 
	
	subplot(1,2,1)
	surfc(RX,RY,(PROB))
	set(gca,'FontSize',16)
	colormap(jet)
	colorbar
	set(gca,'YDir','normal')
%	axis([0 1000 90 200 -1 1])
% 	caxis([-1,1])
	az = -40; 
	el = 35; 
	view(az,el)
	xlabel('Energy [eV]')
	ylabel('Altitude [km]')

	subplot(1,2,2)
	contourf(RX,RY,(PROB))
	set(gca,'FontSize',16)
	colormap(jet)
	colorbar
% 	caxis([-3,2])
%  h = colorbar;
% 	ylabel(h,'LOG10 Ensemble Percentage','FontSize',16);
	xlabel('Energy [eV]')
	ylabel('Altitude [km]')

	fn  = sprintf('./Plots/HvEclick/%04d.jpeg', click); 
	print('-djpeg100',fn); 	
end






