
clear all
close all

fE  = fopen('../planet_X_energy.dat'); 
fH  = fopen('../planet_X_height.dat'); 
fv  = fopen('../planet_X_unit_vel.dat'); 
fd	= fopen('../planet_X_energy_loss.dat'); 
fN	= fopen('../planet_X_Ncoll.dat'); 
RE  = fscanf(fE,'%f',[1,inf]); 
RH  = fscanf(fH,'%f',[1,inf]); 
Rv  = fscanf(fv,'%f',[1,inf]); 
Rd  = fscanf(fd,'%f',[1,inf]); 
RN  = fscanf(fN,'%f',[1,inf]); 
RE  = RE'; 
RH  = RH'/1000;		% convert to km 
Rv  = Rv'; 
Rd  = Rd'; 
RN  = RN'; 

fHE = fopen('../planet_height_vs_energy.dat'); 
fHd = fopen('../planet_height_vs_energy_loss.dat'); 
fHv = fopen('../planet_height_vs_vert_vel.dat'); 
fHN = fopen('../planet_height_vs_Ncoll.dat'); 
dHE = fscanf(fHE,'%d %d %d %f',[4,inf]); 
dHd = fscanf(fHd,'%d %d %d %f',[4,inf]); 
dHv = fscanf(fHv,'%d %d %d %f',[4,inf]); 
dHN = fscanf(fHN,'%d %d %d %f',[4,inf]); 
dHE = dHE'; 
dHd = dHd'; 
dHv = dHv'; 
dHN = dHN'; 

ck  = dHE(:,1); 
Ni  = dHE(:,2); 
Nj  = dHE(:,3); 
PHE = dHE(:,4); 
PHd = dHd(:,4); 
PHv = dHv(:,4); 
PHN = dHN(:,4); 

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
		tot_Hv = 0; 
		tot_HN = 0; 
		for i=1:NI
			for j=1:NJ
				PROB_HE(i,j) = 100*PHE(k); 
				PROB_Hd(i,j) = 100*PHd(k); 
				PROB_Hv(i,j) = 100*PHv(k); 
				PROB_HN(i,j) = 100*PHN(k); 
				tot_HE = tot_HE + PHE(k); 
				tot_Hd = tot_Hd + PHd(k); 
				tot_Hv = tot_Hv + PHv(k); 
				tot_HN = tot_HN + PHN(k); 
				k = k+1; 
			end
		end
		tit = sprintf('%04d',CS*click); 
		suptitle(tit)

		subplot(2,2,1) 
		imagesc(RN,RH,log10(PROB_HN))
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Number of Collisions')
		ylabel('Altitude [km]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_HN*100); 
		title(tit)

		subplot(2,2,2) 
		imagesc(Rv,RH,log10(PROB_Hv))
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Vertical Directional Cosine')
		ylabel('Altitude [km]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_Hv*100); 
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

		fn  = sprintf('./Plots/Height_Command/%04d.jpeg', click); 
		print('-djpeg100',fn); 	

end





