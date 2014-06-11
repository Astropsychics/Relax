
clear all
close all

fx  = fopen('../planet_X_height.dat'); 
fy  = fopen('../planet_X_Ncoll.dat'); 
RX  = fscanf(fx,'%f',[1,inf]); 
RY  = fscanf(fy,'%f',[1,inf]); 
RX  = RX'/1000; 
dx  = RX(2)-RX(1); 
RX  = RX + dx/2; 
RY  = RY'; 

fid = fopen('../planet_height_vs_Ncoll.dat'); 
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
imagesc(RX,RY,PROB)
colormap(hot)
set(gca,'FontSize',16)
set(gca,'YDir','normal')
xlabel('Alitutde [km]')
ylabel('Number of Collisions')

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

figure
plot(RX,AVG_E,'r','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Altitude [km]')
ylabel('Number of Collisions')

figure
plot(RX,STD_E,'g','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Altitude [km]')
ylabel('Standard Deviation')
		
figure
errorbar(RX,AVG_E,STD_L,STD_U,'k*','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Altitude [km]')
ylabel('Average Number of Collisions')
axis([min(RX) max(RX) (min(AVG_E)-max(STD_L)) (max(AVG_E)+max(STD_U))])
print -djpeg100 ./Plots/Averages/AvgN_vs_H.jpeg


