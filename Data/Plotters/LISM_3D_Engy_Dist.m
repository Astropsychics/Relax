
clear all

SMOOTH 				= 0; 
SLICES 				= 0; 
Energy_SLICES = 0; 

AUtoLY = 1.58e-5; 
I0     = 1.1e36; 

f = fopen('../LISM_ENA_3D_Energy_Dist.dat'); 
dat = fscanf(f,'%f %f %f',[3,inf]); 
dat = dat'; 

LD = length(dat(:,1)); 		% Number of rows in dat 
R1 = dat(1,1);						% First distance of dat 

for i=1:LD
	if (dat(i,1) ~= R1)
		LE = i-1;							% Number of angles per energy 
		break
	end
end

NE = LD/LE;								% Number of energies 
NR = NE; 

for i=1:LE
  Engy(i) = dat(i,2); 
end

for i=1:NE
	for j=1:LE
		k = (i-1)*LE + j; 		% Current row count in dat
 		if (dat(k,3) == 0)
			RZ(i,j) = 0; 
			LZ(i,j) = 0; 
		else  
	 		LZ(i,j) = log10(dat(k,3));	% Create matrix, i=Energy, j=Angle, Z(Energy,Angle) = Amp
			RZ(i,j) = dat(k,3); 	
		end
  end
  Distance(i) = dat(k,1); 
end

RZ = RZ'; 
LZ = LZ'; 

FIG_ON = 1; 

if (FIG_ON == 1)
figure
%imagesc(Distance*AUtoLY,Engy,RZ*100)
imagesc(Distance*AUtoLY,Engy,LZ)
colormap(hot)
h = colorbar; 
%ylabel(h,'Percent of Ensemble','FontSize',16);
ylabel(h,'Log_{10} Percent Ensemble','FontSize',16);
set(gca,'YDir','normal'); 
set(gca,'FontSize',16); 
ylabel('ENA Energy [eV]'); 
xlabel('Distance from Star [LY]'); 
print -depsc2 ./Plots/LISM_3D_Energy_Dist.eps
end


if (Energy_SLICES == 1)
	E1 = 500; 
	E2 = 1000; 
	E3 = 1500; 
	E4 = 2000; 
	E5 = 2500; 
	E6 = 3000; 
	E7 = 12000; 

	z1 = 1; 
	while (Engy(z1) <= E1)
		z1 = z1+1; 
	end	
	z2 = 1; 
	while (Engy(z2) <= E2)
		z2 = z2+1; 
	end	
	z3 = 1; 
	while (Engy(z3) <= E3)
		z3 = z3+1; 
	end	
	z4 = 1; 
	while (Engy(z4) <= E4)
		z4 = z4+1; 
	end	
	z5 = 1; 
	while (Engy(z5) <= E5)
		z5 = z5+1; 
	end	
	z6 = 1; 
	while (Engy(z6) <= E6)
		z6 = z6+1; 
	end	
	z7 = 1; 
	while (Engy(z7) <= E7)
		z7 = z7+1; 
	end	

	ES1(1:LE) = 0; 
	for i=1:z1			% loop over all energies in interval
		for j=1:LE		% loop over all distances for energy
			ES1(j) = ES1(j) + RZ(i,j); 	
		end
	end
	ES2(1:LE) = 0; 
	for i=z1:z2			% loop over all energies in interval
		for j=1:LE		% loop over all distances for energy
			ES2(j) = ES2(j) + RZ(i,j); 	
		end
	end
	ES3(1:LE) = 0; 
	for i=z2:z3			% loop over all energies in interval
		for j=1:LE		% loop over all distances for energy
			ES3(j) = ES3(j) + RZ(i,j); 	
		end
	end
	ES4(1:LE) = 0; 
	for i=z3:z4			% loop over all energies in interval
		for j=1:LE		% loop over all distances for energy
			ES4(j) = ES4(j) + RZ(i,j); 	
		end
	end
	ES5(1:LE) = 0; 
	for i=z4:z5			% loop over all energies in interval
		for j=1:LE		% loop over all distances for energy
			ES5(j) = ES5(j) + RZ(i,j); 	
		end
	end
	ES6(1:LE) = 0; 
	for i=z5:z6			% loop over all energies in interval
		for j=1:LE		% loop over all distances for energy
			ES6(j) = ES6(j) + RZ(i,j); 	
		end
	end
	ES7(1:LE) = 0; 
	for i=z6:z7			% loop over all energies in interval
		for j=1:LE		% loop over all distances for energy
			ES7(j) = ES7(j) + RZ(i,j); 	
		end
	end

	AUTOCM = 1.5e13; 

	ES1 = ES1*I0; 
	ES2 = ES2*I0; 
	ES3 = ES3*I0; 
	ES4 = ES4*I0; 
	ES5 = ES5*I0; 
	ES6 = ES6*I0; 
	ES7 = ES7*I0; 

	for i=1:length(ES1)
		ES1(i) = ES1(i)/(Distance(i)*AUTOCM)^2; 
		ES2(i) = ES2(i)/(Distance(i)*AUTOCM)^2; 
		ES3(i) = ES3(i)/(Distance(i)*AUTOCM)^2; 
		ES4(i) = ES4(i)/(Distance(i)*AUTOCM)^2; 
		ES5(i) = ES5(i)/(Distance(i)*AUTOCM)^2; 
		ES6(i) = ES6(i)/(Distance(i)*AUTOCM)^2; 
		ES7(i) = ES7(i)/(Distance(i)*AUTOCM)^2; 
		TOT(i) = ES1(i) + ES2(i) + ES3(i) + ES4(i) + ES5(i) + ES6(i) + ES7(i); 
	end

	DD = Distance*AUtoLY; 

	figure
	semilogy(DD,ES1,'k',DD,ES2,'b',DD,ES3,'g',DD,ES4,'r',DD,ES5,'c', ...
						DD,ES6,'m',DD,ES7,'y','LineWidth',2.5)
	set(gca,'FontSize',16)
	xlabel('Distance from Star [LY]')
	ylabel('ENA Flux [1/cm^2/s]')	
	legend('< 0.5 keV','0.5 keV - 1 keV','1 keV - 1.5 keV','1.5 keV - 2 keV', ...
         '2 keV - 2.5 keV','2.5 keV - 3 keV', '3 keV - 12 keV','Location','Best')

	figure
	semilogy(DD,TOT,'k','LineWidth',2.5)
	set(gca,'FontSize',16)
	xlabel('Distance from Star [LY]')
	ylabel('ENA Flux [1/cm^2/s]')	


end

if (SLICES == 1)

D1 = 1; 
D2 = 5; 
D3 = 10; 
D4 = 15; 
D5 = 25; 

z1 = 1; 
while (Distance(z1)*AUtoLY <= D1)
	z1 = z1+1; 
end
z2 = 1; 
while (Distance(z2)*AUtoLY <= D2)
	z2 = z2+1; 
end
z3 = 1; 
while (Distance(z3)*AUtoLY <= D3)
	z3 = z3+1; 
end
z4 = 1; 
while (Distance(z4)*AUtoLY <= D4)
	z4 = z4+1; 
end
z5 = 1; 
while (Distance(z5)*AUtoLY <= D5)
	z5 = z5+1; 
end

%for i=1:length(Engy)
%	fprintf('Engy(%d): %f\tDist(%d): %f\n', i, Engy(i), i, Distance(i)*AUtoLY); 
%end

if (SMOOTH == 1)

xxx = Engy(1):(Engy(length(Engy))-Engy(1))/100:Engy(length(Engy)); 

for i=1:length(xxx)
	yy1(i) = smoother( Engy, RZ(:,z1), 10, xxx(i) ); 
	yy2(i) = smoother( Engy, RZ(:,z2), 10, xxx(i) ); 
	yy3(i) = smoother( Engy, RZ(:,z3), 10, xxx(i) ); 
	yy4(i) = smoother( Engy, RZ(:,z4), 10, xxx(i) ); 
	yy5(i) = smoother( Engy, RZ(:,z5), 10, xxx(i) ); 
end

figure
plot(xxx,yy1/sum(yy1),'k',xxx,yy2/sum(yy2),'b',xxx,yy3/sum(yy3) ...
,'g',xxx,yy4/sum(yy4),'r','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('ENA Energy [eV]'); 
ylabel('Normalized Frequency'); 
legend('1 LY','10 LY','20 LY','50 LY','Location','Best')
print -depsc2 ./Plots/LISM_Energy_Dist_Smooth_Slices.eps

end

NNN = 1; 

figure
semilogy(Engy,RZ(:,z1)/NNN, 'k', 'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('ENA Energy [eV]'); 
ylabel('Frequency');
title('1 LY From Star') 
print -depsc2 ./Plots/LISM_Energy_Dist_Slices_1LY.eps

figure
semilogy(Engy,RZ(:,z2)/NNN, 'b', 'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('ENA Energy [eV]'); 
ylabel('Frequency');
title('10 LY From Star') 
print -depsc2 ./Plots/LISM_Energy_Dist_Slices_10LY.eps

figure
semilogy(Engy,RZ(:,z3)/NNN, 'g', 'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('ENA Energy [eV]'); 
ylabel('Frequency');
title('20 LY From Star') 
print -depsc2 ./Plots/LISM_Energy_Dist_Slices_20LY.eps

figure
semilogy(Engy,RZ(:,z4)/NNN, 'r', 'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('ENA Energy [eV]'); 
ylabel('Frequency');
title('50 LY From Star') 
print -depsc2 ./Plots/LISM_Energy_Dist_Slices_50LY.eps

figure
semilogy(Engy,RZ(:,z5)/NNN, 'c', 'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('ENA Energy [eV]'); 
ylabel('Frequency');
title('100 LY From Star') 
print -depsc2 ./Plots/LISM_Energy_Dist_Slices_100LY.eps


figure
semilogy(Engy,RZ(:,z1)/NNN, 'k', Engy,RZ(:,z2)/NNN, ...
'b', Engy,RZ(:,z3)/NNN,'g',Engy,RZ(:,z4)/NNN,'r', ...
Engy,RZ(:,z5)/NNN,'c','LineWidth',1.5); 
set(gca,'FontSize',16); 
xlabel('ENA Energy [eV]'); 
ylabel('Normalized Frequency'); 
%axis([10 3000 0 4.5e-3])
legend('1 LY','10 LY','20 LY','50 LY','100 LY','Location','Best')
print -depsc2 ./Plots/LISM_Energy_Dist_Slices.eps

fprintf('SUMS\nz1: %d\tz2: %d\tz3: %d\tz4: %d\tz5: %d\n', sum(RZ(:,z1)), sum(RZ(:,z2)), sum(RZ(:,z3)), sum(RZ(:,z4)), sum(RZ(:,z5))); 

end
