

clear all
close all

f1 = fopen('../planet_all_energy_vs_time.dat'); 
fx = fopen('../planet_X_time.dat'); 
fy = fopen('../planet_X_energy.dat'); 

d1 = fscanf(f1,'%f %f %f',[3,inf]); 
X  = fscanf(fx,'%f',[1,inf]); 
Y  = fscanf(fy,'%f',[1,inf]); 

d1 = d1'; 
X  = X'; 
Y  = Y'; 
dY = Y(2) - Y(1); 
Y  = Y + dY/2; 

Nx = max(d1(:,1)); 
Ny = max(d1(:,2)); 

tot = 0; 
k  = 1; 
for i=1:Nx
	for j=1:Ny
		P(i,j) = d1(k,3); 
		tot    = tot + P(i,j); 
		k = k+1; 
	end
end

P = P'/tot; 

figure
imagesc(X,Y,log10(P))
colormap(hot)
set(gca,'FontSize',16)
set(gca,'Ydir','Normal')
xlabel('Time [Sec]') 
ylabel('Energy [eV]')
colorbar

%slice = [10 15 20 25 40]; 
slice = 1:length(P(:,1)); 

for i=1:length(slice)
	PP    = P(:,slice(i)); 
	YY(i) = sum(PP); 
	WA(i) = 0; 
	for j=1:length(PP)
		WA(i) = WA(i) + PP(j)*Y(j)/YY(i); 	
	end
	figure
	plot(Y,PP/YY(i),'k',WA(i),0,'ro','LineWidth',2.5)
	set(gca,'FontSize',16)
	xlabel('Energy [eV]')
	ylabel('Frequency [1/eV]')
	tit = sprintf('%f sec slice', X(slice(i))); 
	title(tit)
end

ALL_PLOT_ON = 0; 
if (ALL_PLOT_ON == 1)
figure
plot(Y,P(:,slice(1))/YY(1),'k', ...
     Y,P(:,slice(2))/YY(2),'b', ...
     Y,P(:,slice(3))/YY(3),'g', ...
     Y,P(:,slice(4))/YY(4),'r', ...
     Y,P(:,slice(5))/YY(5),'c', ...
		'LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Energy [eV]')
ylabel('Frequency [1/eV]')
end

for i=1:length(P(:,1))
	PP 		= P(:,i); 
	YY(i) = sum(PP); 
	WA(i) = 0; 
	for j=1:length(PP)
		WA(i) = WA(i) + PP(j)*Y(j)/YY(i); 
	end
end

figure
plot(X,WA,'r.',X,WA,'k','LineWidth',2.5)
set(gca,'FontSize',16)
ylabel('Weighted Average Energy [eV]')
xlabel('Time [sec]')


