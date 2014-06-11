
clear all
close all

fE  = fopen('../lism_X_energy.dat'); 
fX  = fopen('../lism_X_xyz.dat'); 
fR  = fopen('../lism_X_r.dat'); 
fT  = fopen('../lism_X_time.dat'); 

dE  = fscanf(fE,'%f',[1,inf]); 
dX  = fscanf(fX,'%f',[1,inf]); 
dR  = fscanf(fR,'%f',[1,inf]); 
dT  = fscanf(fT,'%f',[1,inf]); 

dE  = dE'; 
dX  = dX'; 
dR  = dR'; 
dT  = dT'; 

sec_to_year = 1/31536000;
AU_to_PC    = 4.84813681e-6;   
dT  = dT*sec_to_year;
dX  = dX*AU_to_PC; 
dR  = dR*AU_to_PC; 

fET = fopen('../lism_energy_vs_time.dat'); 
fXT = fopen('../lism_X_vs_time.dat'); 
fYT = fopen('../lism_Y_vs_time.dat'); 
fZT = fopen('../lism_Z_vs_time.dat'); 
fRT = fopen('../lism_R_vs_time.dat'); 
dET = fscanf(fET,'%d %d %d %f',[4,inf]); 
dXT = fscanf(fXT,'%d %d %d %f',[4,inf]); 
dYT = fscanf(fYT,'%d %d %d %f',[4,inf]); 
dZT = fscanf(fZT,'%d %d %d %f',[4,inf]); 
dRT = fscanf(fRT,'%d %d %d %f',[4,inf]); 

dET = dET'; 
dXT = dXT'; 
dYT = dYT'; 
dZT = dZT'; 
dRT = dRT'; 

ck  = dRT(:,1); 
Ni  = dRT(:,2); 
Nj  = dRT(:,3); 
PET = dET(:,4); 
PXT = dXT(:,4); 
PYT = dYT(:,4); 
PZT = dZT(:,4); 
PRT = dRT(:,4); 

N   = length(PET); 

NI  = max(Ni); 
NJ  = max(Nj); 
NC  = N/(NI*NJ); 
CS  = max(ck)/NC; 

k = 1; 
figure('units','normalized','outerposition',[0 0 1 1])
for click=1:NC
		tot_XT = 0; 
		tot_YT = 0; 
		tot_ZT = 0; 
		tot_RT = 0; 
		for i=1:NI
			for j=1:NJ
				PROB_XT(i,j) = 100*PXT(k); 
				PROB_YT(i,j) = 100*PYT(k); 
				PROB_ZT(i,j) = 100*PZT(k); 
				PROB_RT(i,j) = 100*PRT(k); 
				tot_XT = tot_XT + PXT(k); 
				tot_YT = tot_YT + PYT(k); 
				tot_ZT = tot_ZT + PZT(k); 
				tot_RT = tot_RT + PRT(k); 
				k = k+1; 
			end
		end
		tit = sprintf('%04d',CS*click); 
		suptitle(tit)

		subplot(2,2,1) 
		contourf(dT,dX,PROB_XT')
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Time [days]')
		ylabel('X [PC]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_XT*100); 
		title(tit)

		subplot(2,2,2) 
		contourf(dT,dX,PROB_YT')
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Time [days]')
		ylabel('Y [PC]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_YT*100); 
		title(tit)

		subplot(2,2,3) 
		contourf(dT,dX,PROB_ZT')
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Time [days]')
		ylabel('Z [PC]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_ZT*100); 
		title(tit)

		subplot(2,2,4) 
		contourf(dT,dR,PROB_RT')
		colormap(hot)
		set(gca,'FontSize',12)
		set(gca,'YDir','normal')
		xlabel('Time [days]')
		ylabel('R [PC]')
		tit = sprintf('%5.2f Percent Ensemble in Frame',tot_RT*100); 
		title(tit)
	  h = colorbar('EastOutside');
		set(h,'Position',[0.92 0.11 0.01 0.82])
 	 	ylabel(h,'Ensemble Percentage','FontSize',10);
 	 	caxis([-3,2])

		fn  = sprintf('./Plots/LISM_PosTime_Command/%04d.jpeg', click); 
		print('-djpeg100',fn); 	

end





