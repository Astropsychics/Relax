
clear all
close all

fE  = fopen('../planet_X_energy.dat'); 
fH  = fopen('../planet_X_height.dat'); 
fT  = fopen('../planet_X_time.dat'); 
fd	= fopen('../planet_X_energy_loss.dat'); 
RE  = fscanf(fE,'%f',[1,inf]); 
RH  = fscanf(fH,'%f',[1,inf]); 
RT  = fscanf(fT,'%f',[1,inf]); 
Rd  = fscanf(fd,'%f',[1,inf]); 
RE  = RE'; 
RH  = RH'/1000;		% convert to km 
RT  = RT'; 
Rd  = Rd'; 

fHE = fopen('../planet_height_vs_energy.dat'); 
fHd = fopen('../planet_height_vs_energy_loss.dat'); 
fHT = fopen('../planet_height_vs_time.dat'); 
fET = fopen('../planet_energy_vs_time.dat'); 
dHE = fscanf(fHE,'%d %d %d %f',[4,inf]); 
dHd = fscanf(fHd,'%d %d %d %f',[4,inf]); 
dHT = fscanf(fHT,'%d %d %d %f',[4,inf]); 
dET = fscanf(fET,'%d %d %d %f',[4,inf]); 
dHE = dHE'; 
dHd = dHd'; 
dHT = dHT'; 
dET = dET'; 

ck  = dHE(:,1); 
Ni  = dHE(:,2); 
Nj  = dHE(:,3); 
PHE = dHE(:,4); 
PHd = dHd(:,4); 
PHT = dHT(:,4); 
PET = dET(:,4); 

N   = length(PHE); 

NI  = max(Ni); 
NJ  = max(Nj); 
NC  = N/(NI*NJ); 
CS  = max(ck)/NC;  

k = 1; 
figure('units','normalized','outerposition',[0 0 1 1])
for click=1:NC
		tot_HE = 0; 
		tot_Hd = 0; 
		tot_HT = 0; 
		tot_ET = 0; 
		for i=1:NI
			for j=1:NJ
				PROB_HE(i,j) = 100*PHE(k); 
				PROB_Hd(i,j) = 100*PHd(k); 
				PROB_HT(i,j) = 100*PHT(k); 
				PROB_ET(i,j) = 100*PET(k); 
				tot_HE = tot_HE + PHE(k); 
				tot_Hd = tot_Hd + PHd(k); 
				tot_HT = tot_HT + PHT(k); 
				tot_ET = tot_ET + PET(k); 
				k = k+1; 
			end
		end
		tit = sprintf('%04d',CS*click); 
		suptitle(tit)

		subplot(2,2,1) 
		imagesc(RT,RE,log10(PROB_ET))
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Time [sec]')
		ylabel('Energy [eV]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_ET*100); 
		title(tit)
%	  h = colorbar('EastOutside');
% 	 	ylabel(h,'LOG10 Ensemble Percentage','FontSize',10);
% 	 	caxis([-3,2])

		subplot(2,2,2) 
		imagesc(RT,RH,log10(PROB_HT))
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Time [sec]')
		ylabel('Altitude [km]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_HT*100); 
		title(tit)

		subplot(2,2,3) 
		imagesc(RE,RH,log10(PROB_HE))
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Energy [eV]')
		ylabel('Altitude [km]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_HE*100); 
		title(tit)

		subplot(2,2,4) 
		imagesc(Rd,RH,log10(PROB_Hd))
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Energy Loss [eV]')
		ylabel('Altitude [km]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_Hd*100); 
		title(tit)

	  h = colorbar('EastOutside');
		set(h,'Position',[0.92 0.11 0.01 0.82])
 	 	ylabel(h,'LOG10 Ensemble Percentage','FontSize',10);
 	 	caxis([-3,2])

		fn  = sprintf('./Plots/Command/%04d.jpeg', click); 
		print('-djpeg100',fn); 	

end





