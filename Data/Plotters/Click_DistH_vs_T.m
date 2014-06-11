
clear all
close all

fx  = fopen('../planet_X_time.dat'); 
fy  = fopen('../planet_X_height.dat'); 
RX  = fscanf(fx,'%f',[1,inf]); 
RY  = fscanf(fy,'%f',[1,inf]); 
RX  = RX'; 
dx  = RX(2)-RX(1); 
RX  = RX + dx/2; 
RY  = RY'/1000; 

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
NC  = N/(NI*NJ); 
CS  = max(ck)/NC; 

PROB(1:NI,1:NJ) = 0; 

k   = 1; 
tot = 0; 

CLICK_COUNT(1:NJ) = 0; 

for click=1:NC
	for i=1:NI
		for j=1:NJ
			P_NOW(i,j) = P(k); 
			tot	= tot + P(k);
			k = k+1; 
		end
	end
	PROB(:,:) = PROB(:,:) + P_NOW(:,:); 
	for j=1:NJ
		if ( sum(P_NOW(:,j)) > 0 )
			CLICK_COUNT(j) = CLICK_COUNT(j) + 1; 
		end
	end
end

figure
imagesc(RX,RY,PROB/tot)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','normal')
xlabel('Time [sec]')
ylabel('Altitude [km]')

dE = RY(2) - RY(1); 
figure

for i=1:NI

	T(i)   = sum(PROB(:,i)); 
	CCC    = T(i)/CLICK_COUNT(i); 
	Y(i,:) = PROB(:,i);
%	Y(i,:) = PROB(:,i)/CLICK_COUNT(i);
	TTT(i) = sum(Y(i,:))/tot;  
%	Y(i,:) = PROB(:,i)/NC; 

	plot(RY,Y(i,:),'g*',RY,Y(i,:),'k','LineWidth',2.5)
	set(gca,'FontSize',16)
	xlabel('Altitude [km]')
	ylabel('Ensemble Percentage')
	axis([70 200 0 100])
	tit = sprintf('%f < t < %f sec\n%f Percent Ensemble in Frame\n%d Clicks Contributing', RX(i)-dx/2, RX(i)+dx/2,TTT(i),CLICK_COUNT(i)); 
	title(tit); 
	fn  = sprintf('./Plots/Hdist_T/%03d_Hdist_T.jpeg',i); 
	print('-djpeg100',fn) 
	
end



