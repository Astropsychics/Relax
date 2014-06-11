
clear all
close all

fx = fopen('../lism_X_time.dat'); 
fy = fopen('../lism_X_energy.dat'); 
RX  = fscanf(fx,'%f',[1,inf]); 
RY  = fscanf(fy,'%f',[1,inf]); 
RX  = RX'; 
dx  = RX(2)-RX(1); 
RY  = RY'; 
dy  = RY(2)-RY(1); 

fid = fopen('../lism_all_energy_vs_time.dat'); 
dat = fscanf(fid,'%d %d %f',[3,inf]); 
dat = dat'; 

Ni  = dat(:,1); 
Nj  = dat(:,2); 
P   = dat(:,3); 

N   = length(P); 

NI  = max(Ni); 
NJ  = max(Nj); 
NC  = N/(NI*NJ); 

PROB(1:NI,1:NJ) = 0; 

k   = 1; 
tot = 0; 

for click=1:NC
	for i=1:NI
		for j=1:NJ
			PROB(i,j) = PROB(i,j) + P(k); 
			tot = tot + P(k); 
			k = k+1; 
		end
	end
end

TOTAL = sum(sum(PROB(:,:))*dx*dy); 
PROB  = PROB/TOTAL; 

year = pi*1e7; 
RX   = RX/year; 
dx   = RX(2) - RX(1); 
dy   = RY(2) - RY(1); 

C1 = 1/(sum(sum(PROB(:,:)))*dx*dy); 

PROB = PROB*C1; 

figure
imagesc(RX,RY,PROB)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','normal')
xlabel('Time [years]'); 
ylabel('Energy [eV]'); 
axis([5000 30000 0 12000])
h = colorbar;
ylabel(h,'Probability Density [1/year/eV]','FontSize',16);
print -depsc2 ./Plots/LISM_All_EvT_PD.eps

dE = RY(2) - RY(1); 
for i=1:NI
	BOT = 0;
	TOP = 0;  
	for j=1:NJ
		BOT  = BOT + PROB(j,i); 
		TOP  = TOP + RY(j)*PROB(j,i); 
	end
	AVG_E(i) = TOP/BOT; 

	UT = 0; 
	UB = 0; 
	LT = 0; 
	LB = 0; 

	kU = 1; 
	kL = 1; 
	for j=1:NJ
		if ( RY(j) >= AVG_E(i) )
			% upper std	
			UT = UT + PROB(j,i)*(RY(j)-AVG_E(i))^2; 
			UB = UB + PROB(j,i); 	
			kU = kU + 1; 
		else
			% lower std	
			LT = LT + PROB(j,i)*(RY(j)-AVG_E(i))^2; 
			LB = LB + PROB(j,i); 	
			kL = kL + 1; 
		end
	end

	STD_U(i) = sqrt(UT/UB); 
	STD_L(i) = sqrt(LT/LB); 

	TOP = 0; 
	for j=1:NJ
		TOP = TOP + PROB(j,i)*(RY(j)-AVG_E(i))^2; 
	end
	STD_E(i) = sqrt(TOP/BOT); 	
end

MEAN_STD_ON = 0; 

if (MEAN_STD_ON == 1)
figure
plot(RX,AVG_E,'r','LineWidth',2.5)
set(gca,'FontSize',16)
title(tit); 
xlabel('X')
ylabel('Average Y')

figure
plot(RX,STD_E,'g','LineWidth',2.5)
set(gca,'FontSize',16)
title(tit); 
xlabel('X')
ylabel('Standard Deviation Y')

figure
errorbar(RX,AVG_E,STD_L,STD_U,'k*','LineWidth',2.5)
set(gca,'FontSize',16)
title(tit); 
xlabel(xlb); 
ylabel('Average Y'); 
%axis([min(RX) max(RX) (min(AVG_E)-max(STD_L)) (max(AVG_E)+max(STD_U))])
fnam = sprintf('./Plots/Averages/%s.jpeg', NM); 
print('-djpeg100',fnam); 
%print -djpeg100 ./Plots/Averages/AvgE_vs_T.jpeg
end

