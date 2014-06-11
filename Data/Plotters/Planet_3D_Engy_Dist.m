
clear all

SLICES 				= 0; 
SMOOTH 				= 0; 

%I0     				= 1e6; % ENA's/sec incident on an area of (R/4)^2 R=mars Radius 
%I0     				= 33.3; % ENA's/sec incident on an area of (R/4)^2 R=mars Radius 
%I0     				= 7.2e21; % ENA's/sec incident on an area of (R/4)^2 R=mars Radius 
												% Use incident flux from Kallio et al. 1997
												% J(sza=0) = 10^6 ENAs/cm^2/sec

f = fopen('../planet_ENA_3D_Energy_Dist.dat'); 
dat = fscanf(f,'%f %f %f',[3,inf]);  % [ height [m] | energy [eV] | # particles ] 
dat = dat'; 

LD = length(dat(:,1)); 		% Number of rows in dat 
R1 = dat(1,1);						% First height of dat 

for i=1:LD
	if (dat(i,1) ~= R1)
		LE = i-1;							% Number of energies per height
		break
	end
end

NE = LD/LE;								% Number of heights
NR = NE; 

for i=1:LE
  Engy(i) = dat(i,2); 
end

Q = 0; 

for i=1:NE
	for j=1:LE
		k = (i-1)*LE + j; 
		Q = Q + dat(k,3); 
	end
end

I0 = 1/Q; 

for i=1:NE								% loop over heights
	for j=1:LE							% loop over energies
		k = (i-1)*LE + j; 		% Current row count in dat
 		if (dat(k,3) == 0)
			RZ(i,j) = 0; 
			LZ(i,j) = 0; 
		else  
	 		LZ(i,j) = log10(I0*dat(k,3));	% Create matrix, i=height, j=Energy, Z(i,j) = count 
			RZ(i,j) = I0*dat(k,3); 	
		end
  end
  Height(i) = dat(k,1); 
end

AVG_ON = 1; 
if (AVG_ON == 1) 
	for i=1:NE
		for j=1:LE
			ZZ(i,j) = RZ(i,j)*Engy(j); 
		end	
		NN_E(i) = sum(RZ(i,:)); 
	end
	for i=1:NE
		AVG_E(i) = sum(ZZ(i,:))/NN_E(i); 
		STD_E(i) = std(ZZ(i,:))/NN_E(i); 		
	end
	figure
	semilogx(AVG_E,Height/1.0D3,'g','LineWidth',2.5)
	set(gca,'FontSize',16)
	xlabel('Average ENA Energy [eV]')
	ylabel('Altitude [km]')

	figure
	semilogx(STD_E,Height/1.0D3,'r','LineWidth',2.5)
	set(gca,'FontSize',16)
	xlabel('Standard Deviation from Average [eV]')
	ylabel('Altitude [km]')
end

FIG_ON = 1; 

if (FIG_ON == 1)
figure
%imagesc(Engy,Height/1.0D3,RZ)
imagesc(Engy/1.0D3,Height/1.0D3,LZ)
colormap(hot)
h = colorbar; 
ylabel(h,'Log_{10}   cm^{-2} sec^{-1}','FontSize',16);
set(gca,'YDir','normal'); 
set(gca,'FontSize',16); 
xlabel('ENA Energy [keV/amu]'); 
ylabel('Altitude [km]'); 
print -depsc2 ./Plots/Planet_3D_Energy_Dist.eps

figure
imagesc(Engy,Height/1.0D3,100*RZ)
colormap(hot)
h = colorbar; 
ylabel(h,'cm^{-2} sec^{-1}','FontSize',16);
set(gca,'YDir','normal'); 
set(gca,'FontSize',16); 
xlabel('ENA Energy [keV/amu]'); 
ylabel('Altitude [km]'); 
print -depsc2 ./Plots/Planet_3D_Energy_Dist_Percent.eps
end

if (SLICES == 1)

	% Heights to slice up [km]
	D0 = 80; 
	D1 = 85; 
	D2 = 90; 
	D3 = 100; 
	D4 = 125; 
	D5 = 150; 

	MtoKM = 1/1000; 

	z0 = 1; 
	while (Height(z0)*MtoKM <= D0)
		z0 = z0+1; 
	end
	z1 = 1; 
	while (Height(z1)*MtoKM <= D1)
		z1 = z1+1; 
	end
	z2 = 1; 
	while (Height(z2)*MtoKM <= D2)
		z2 = z2+1; 
	end
	z3 = 1; 
	while (Height(z3)*MtoKM <= D3)
		z3 = z3+1; 
	end
	z4 = 1; 
	while (Height(z4)*MtoKM <= D4)
		z4 = z4+1; 
	end
	z5 = 1; 
	while (Height(z5)*MtoKM <= D5)
		z5 = z5+1; 
	end

	if (SMOOTH == 1)

		xxx = Engy(1):(Engy(length(Engy))-Engy(1))/100:Engy(length(Engy)); 

		for i=1:length(xxx)
			yy0(i) = smoother( Engy, RZ(z0,:), 10, xxx(i) ); 
			yy1(i) = smoother( Engy, RZ(z1,:), 10, xxx(i) ); 
			yy2(i) = smoother( Engy, RZ(z2,:), 10, xxx(i) ); 
			yy3(i) = smoother( Engy, RZ(z3,:), 10, xxx(i) ); 
			yy4(i) = smoother( Engy, RZ(z4,:), 10, xxx(i) ); 
			yy5(i) = smoother( Engy, RZ(z5,:), 10, xxx(i) ); 
		end

		figure
		semilogy(xxx,yy0/sum(yy0),'y',xxx,yy1/sum(yy1),'k',xxx,yy2/sum(yy2), ...
		'b',xxx,yy3/sum(yy3),'g',xxx,yy4/sum(yy4),'r',xxx,yy5/sum(yy5),'c','LineWidth',2.5); 
		set(gca,'FontSize',16); 
		xlabel('ENA Energy [eV]'); 
		ylabel('Normalized Frequency'); 
		legend('80 km','85 km','90 km','100 km','125 km','150 km','Location','Best')
		axis([min(xxx) max(xxx) 1e-4 1])
		print -depsc2 ./Plots/Planet_Energy_Dist_Smooth_ALL.eps

	end % SMOOTH

	NNN = 1; 

	figure
	semilogy(Engy,RZ(z0,:)/NNN, 'y', 'LineWidth',2.5); 
	set(gca,'FontSize',16); 
	xlabel('ENA Energy [eV]'); 
	ylabel('Frequency');
	tit = sprintf('%d km', D0); 
	title(tit) 
	print -depsc2 ./Plots/Planet_Energy_Dist_80km.eps

	figure
	semilogy(Engy,RZ(z1,:)/NNN, 'k', 'LineWidth',2.5); 
	set(gca,'FontSize',16); 
	xlabel('ENA Energy [eV]'); 
	ylabel('Frequency');
	tit = sprintf('%d km', D1); 
	title(tit) 
	print -depsc2 ./Plots/Planet_Energy_Dist_85km.eps

	figure
	semilogy(Engy,RZ(z2,:)/NNN, 'b', 'LineWidth',2.5); 
	set(gca,'FontSize',16); 
	xlabel('ENA Energy [eV]'); 
	ylabel('Frequency');
	tit = sprintf('%d km', D2); 
	title(tit) 
	print -depsc2 ./Plots/Planet_Energy_Dist_90km.eps

	figure
	semilogy(Engy,RZ(z3,:)/NNN, 'g', 'LineWidth',2.5); 
	set(gca,'FontSize',16); 
	xlabel('ENA Energy [eV]'); 
	ylabel('Frequency');
	tit = sprintf('%d km', D3); 
	title(tit) 
	print -depsc2 ./Plots/Planet_Energy_Dist_100km.eps

	figure
	semilogy(Engy,RZ(z4,:)/NNN, 'r', 'LineWidth',2.5); 
	set(gca,'FontSize',16); 
	xlabel('ENA Energy [eV]'); 
	ylabel('Frequency');
	tit = sprintf('%d km', D4); 
	title(tit) 
	print -depsc2 ./Plots/Planet_Energy_Dist_125km.eps

	figure
	semilogy(Engy,RZ(z5,:)/NNN, 'c', 'LineWidth',2.5); 
	set(gca,'FontSize',16); 
	xlabel('ENA Energy [eV]'); 
	ylabel('Frequency');
	tit = sprintf('%d km', D5); 
	title(tit) 
	print -depsc2 ./Plots/Planet_Energy_Dist_150km.eps

	figure
	semilogy(Engy,RZ(z0,:)/NNN,'y',Engy,RZ(z1,:)/NNN, 'k', Engy,RZ(z2,:)/NNN, ...
	'b', Engy,RZ(z3,:)/NNN,'g',Engy,RZ(z4,:)/NNN,'r', ...
	Engy,RZ(z5,:)/NNN,'c','LineWidth',1.5); 
	set(gca,'FontSize',16); 
	xlabel('ENA Energy [eV]'); 
	ylabel('Normalized Frequency'); 
	%axis([10 3000 0 4.5e-3])
	legend('80 km','85 km','90 km','100 km','125 km','150 km','Location','Best')
	print -depsc2 ./Plots/Planet_Energy_Dist_All.eps

  PRINT_DATA = 1;
  if (PRINT_DATA == 1)
    fo0 = fopen('./Plots/ENA_80km.dat','w');
    fo1 = fopen('./Plots/ENA_85km.dat','w');
    fo2 = fopen('./Plots/ENA_90km.dat','w');
    fo3 = fopen('./Plots/ENA_100km.dat','w');
    fo4 = fopen('./Plots/ENA_125km.dat','w');
    fo5 = fopen('./Plots/ENA_150km.dat','w');
    for i=1:length(Engy)
      fprintf(fo0,'%e %e\n', Engy(i)/1000, RZ(z0,i)); 
      fprintf(fo1,'%e %e\n', Engy(i)/1000, RZ(z1,i));
      fprintf(fo2,'%e %e\n', Engy(i)/1000, RZ(z2,i));
      fprintf(fo3,'%e %e\n', Engy(i)/1000, RZ(z3,i));
      fprintf(fo4,'%e %e\n', Engy(i)/1000, RZ(z4,i));
      fprintf(fo5,'%e %e\n', Engy(i)/1000, RZ(z5,i));
    end
  end
end % SLICES

