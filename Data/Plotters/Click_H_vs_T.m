
clear all
close all

fx  = fopen('../planet_X_time.dat'); 
fy  = fopen('../planet_X_height.dat'); 
RX  = fscanf(fx,'%f',[1,inf]); 
RY  = fscanf(fy,'%f',[1,inf]); 
RX  = RX'; 
RY  = RY'; 
RY  = RY/1000; 

fid = fopen('../planet_height_vs_time.dat'); 
dat = fscanf(fid,'%d %d %d %f',[4,inf]); 
dat = dat'; 


ck  = dat(:,1); 
Ni  = dat(:,2); 
Nj  = dat(:,3); 
P   = dat(:,4); 

N   = length(P); 

NI  = max(Ni); 
NJ  = max(Nj); 
NC  = max(ck); 

k = 1; 
figure
for click=1:NC

	if (mod(click,10) == 0)
		tot = 0; 
		for i=1:NI
			for j=1:NJ
				PROB(i,j) = 100*P(k); 
				tot = tot + P(k); 
				k = k+1; 
			end
		end
		imagesc(RX,RY,log10(PROB))
		colormap(hot)
	  h = colorbar;
 	 	ylabel(h,'LOG10 Ensemble Percentage','FontSize',16);
 	 	caxis([-3,2])
		set(gca,'FontSize',16)
		set(gca,'YDir','normal')
		tit = sprintf('%04d CLICK\t%5.2f Percent Ensemble in Frame',click,tot*100); 
		title(tit)
		xlabel('Time [sec]')
		ylabel('Altitude [km]')
		fn  = sprintf('./Plots/HvTclick/%04d.jpeg', click); 
		print('-djpeg100',fn); 	
	end

end





