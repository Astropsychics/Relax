
clear all
close all

fx  = fopen('../lism_X_time.dat'); 
fy  = fopen('../lism_X_energy.dat'); 
RX  = fscanf(fx,'%f',[1,inf]); 
RY  = fscanf(fy,'%f',[1,inf]); 
RX  = RX'; 
RY  = RY'; 

fid = fopen('../lism_energy_vs_time.dat'); 
dat = fscanf(fid,'%d %d %d %f',[4,inf]); 
dat = dat'; 

year = pi*1e7;
RX   = RX/year;  

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
figure
for click=1:NC
	tot = 0; 
	for i=1:NI
		for j=1:NJ
			PROB(i,j) = 100*P(k); 
			tot = tot + P(k); 
			k = k+1; 
		end
	end

	PROB = PROB';

	imagesc(RX,RY,log10(PROB))
	colormap(hot)
  h = colorbar;
 	ylabel(h,'LOG10 Ensemble Percentage','FontSize',16);
 	caxis([-3,2])
	set(gca,'FontSize',16)
	set(gca,'YDir','normal')
	tit = sprintf('%04d CLICK\t%5.2f Percent Ensemble in Frame',click,tot*100); 
	title(tit)
	xlabel('Time [year]')
	ylabel('Energy [eV]')
	fn  = sprintf('./Plots/EvTclick/%04d.jpeg', click); 
	print('-djpeg100',fn); 	
end

