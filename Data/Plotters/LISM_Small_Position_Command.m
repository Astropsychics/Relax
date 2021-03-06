
clear all
close all

fX  = fopen('../lism_small_X_xyz.dat'); 
fR  = fopen('../lism_X_energy.dat'); 
%fR  = fopen('../lism_small_X_r.dat'); 
fT  = fopen('../lism_X_time.dat'); 

dX  = fscanf(fX,'%f',[1,inf]); 
dR  = fscanf(fR,'%f',[1,inf]); 
dT  = fscanf(fT,'%f',[1,inf]); 

dX  = dX'; 
dR  = dR'; 
dT  = dT'; 

sec_to_year = 1/31536000;
AU_to_PC    = 4.84813681e-6;   
dT  = dT*sec_to_year;
%dX  = dX*AU_to_PC; 
%dR  = dR*AU_to_PC; 

ddX = dX(2) - dX(1); 
ddR = dR(2) - dR(1); 
ddT = dT(2) - dT(1); 

fXY = fopen('../lism_small_X_vs_Y.dat'); 
fYZ = fopen('../lism_small_Y_vs_Z.dat'); 
fXZ = fopen('../lism_small_X_vs_Z.dat'); 
fRT = fopen('../lism_energy_vs_time.dat'); 
%fRT = fopen('../lism_small_R_vs_time.dat'); 
dXY = fscanf(fXY,'%d %d %d %f',[4,inf]); 
dYZ = fscanf(fYZ,'%d %d %d %f',[4,inf]); 
dXZ = fscanf(fXZ,'%d %d %d %f',[4,inf]); 
dRT = fscanf(fRT,'%d %d %d %f',[4,inf]); 

dXY = dXY'; 
dYZ = dYZ'; 
dXZ = dXZ'; 
dRT = dRT'; 

ck  = dRT(:,1); 
Ni  = dRT(:,2); 
Nj  = dRT(:,3); 
PXY = dXY(:,4); 
PYZ = dYZ(:,4); 
PXZ = dXZ(:,4); 
PRT = dRT(:,4); 

N   = length(PRT); 

NI  = max(Ni); 
NJ  = max(Nj); 
NC  = N/(NI*NJ); 
CS  = max(ck)/NC; 

NNN = 200; 
for i=1:NNN
	theta(i)   = (i-1)*2*pi/(NNN-1); 	
	EARTH_X(i) = cos(theta(i)); 
	EARTH_Y(i) = sin(theta(i)); 
end

k = 1; 
figure('units','normalized','outerposition',[0 0 1 1])
for click=1:NC
		tot_1 = 0; 
		tot_2 = 0; 
		tot_3 = 0; 
		tot_4 = 0; 
		for i=1:NI
			for j=1:NJ
				PROB_1(i,j) = 100*PXY(k); 
				PROB_2(i,j) = 100*PXZ(k); 
				PROB_3(i,j) = 100*PYZ(k); 
				PROB_4(i,j) = 100*PRT(k); 
				tot_1 = tot_1 + PXY(k); 
				tot_2 = tot_2 + PXZ(k); 
				tot_3 = tot_3 + PYZ(k); 
				tot_4 = tot_4 + PRT(k); 
				k = k+1; 
			end
		end

		PROB_1 = PROB_1/(tot_1*ddX*ddX); 
		PROB_2 = PROB_2/(tot_2*ddX*ddX); 
		PROB_3 = PROB_3/(tot_3*ddX*ddX); 
		PROB_4 = PROB_4/(tot_4*ddR*ddT); 

		CB_E = 10; 

		tit = sprintf('%04d',CS*click); 
		suptitle(tit)

		subplot(2,2,1) 
		imagesc(dX,dX,PROB_1')
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('X [AU]')
		ylabel('Y [AU]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_1*100); 
		title(tit)
	  h = colorbar('EastOutside');
 	 	ylabel(h,'Probability Density [1/AU^2]','FontSize',10);
% 	 	caxis([0,CB_E])

		subplot(2,2,2) 
		imagesc(dX,dX,PROB_2')
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('X [AU]')
		ylabel('Z [AU]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_2*100); 
		title(tit)
	  h = colorbar('EastOutside');
 	 	ylabel(h,'Probability Density [1/AU^2]','FontSize',10);
% 	 	caxis([0,CB_E])

		subplot(2,2,3) 
		imagesc(dX,dX,PROB_3')
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Y [AU]')
		ylabel('Z [AU]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_3*100); 
		title(tit)
	  h = colorbar('EastOutside');
 	 	ylabel(h,'Probability Density [1/AU^2]','FontSize',10);
% 	 	caxis([0,CB_E])

		subplot(2,2,4) 
		imagesc(dT,dR,PROB_4')
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Time [years]')
		ylabel('Energy [eV]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_4*100); 
		title(tit)
	  h = colorbar('EastOutside');
 	 	ylabel(h,'Probability Density [1/years/eV]','FontSize',10);
% 	 	caxis([0,1e-3])

%	  h = colorbar('EastOutside');
%		set(h,'Position',[0.92 0.11 0.01 0.82])
% 	 	ylabel(h,'Ensemble Percentage','FontSize',10);
% 	 	caxis([-3,2])

		fn  = sprintf('./Plots/LISM_Small_Position_Command/%04d.jpeg', click); 
		print('-djpeg100',fn); 	

end





