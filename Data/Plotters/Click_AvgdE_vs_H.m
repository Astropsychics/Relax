
clear all
close all

fy  = fopen('../planet_X_height.dat'); 
fx  = fopen('../planet_X_energy_loss.dat'); 
RX  = fscanf(fx,'%f',[1,inf]); 
RY  = fscanf(fy,'%f',[1,inf]); 
RX  = RX'; 
dx  = RX(2)-RX(1); 
RX  = RX + dx/2; 
RY  = RY'/1000; 

fid = fopen('../planet_height_vs_energy_loss.dat'); 
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

for click=1:NC
	for i=1:NI
		for j=1:NJ
			PROB(i,j) = PROB(i,j) + P(k); 
			tot = tot + P(k); 
			k = k+1; 
		end
	end
end

PROB(:,:) = PROB(:,:)/NC; 

figure
imagesc(RY,RX,PROB')
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','normal')
ylabel('Energy Loss [eV]')
xlabel('Altitude [km]')

dZ = RY(2) - RY(1); 

for i=1:NI
	BOT = 0;
	TOP = 0;  
	for j=1:NJ
		BOT  = BOT + PROB(i,j); 
		TOP  = TOP + RX(j)*PROB(i,j); 
	end
	AVG_E(i) = TOP/BOT; 

	UT = 0; 
	UB = 0; 
	LT = 0; 
	LB = 0; 

	kU = 1; 
	kL = 1; 
	for j=1:NJ
		if ( RX(j) >= AVG_E(i) )
			% upper std	
			UT = UT + PROB(i,j)*(RX(j)-AVG_E(i))^2; 
			UB = UB + PROB(i,j); 	
			kU = kU + 1; 
		else
			% lower std	
			LT = LT + PROB(i,j)*(RX(j)-AVG_E(i))^2; 
			LB = LB + PROB(i,j); 	
			kL = kL + 1; 
		end
	end

	STD_U(i) = sqrt(UT/UB); 
	STD_L(i) = sqrt(LT/LB); 

	TOP = 0; 
	for j=1:NJ
		TOP = TOP + PROB(i,j)*(RX(j)-AVG_E(i))^2; 
	end
	STD_E(i) = sqrt(TOP/BOT); 	
end

figure
plot(RY,AVG_E,'r','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Altitude [km]')
ylabel('Average Energy Loss [eV]')

figure
plot(RY,STD_E,'g','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Altitude [km]')
ylabel('Standard Deviation [eV]')
		
figure
errorbar(RY,AVG_E,STD_L,STD_U,'k*','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Altitude [km]')
ylabel('Average Energy Loss [eV]')
axis([min(RY) max(RY) (min(AVG_E)-max(STD_L)) (max(AVG_E)+max(STD_U))])
print -djpeg100 ./Plots/Averages/AvgdE_vs_H.jpeg


