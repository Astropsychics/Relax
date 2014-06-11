
clear all
close all

f_heo    = fopen('../HeO_Average_Scattering_Angles.dat'); 
f_hehe   = fopen('../HeHe_Average_Scattering_Angles.dat'); 
f_heh    = fopen('../HeH_Average_Scattering_Angles.dat'); 
f_hh     = fopen('../HH_Average_Scattering_Angles.dat'); 
f_ho     = fopen('../HO_Average_Scattering_Angles.dat'); 
f_har    = fopen('../HAr_Average_Scattering_Angles.dat'); 
f_hh2    = fopen('../HH2_Average_Scattering_Angles.dat'); 
f_hn2    = fopen('../HN2_Average_Scattering_Angles.dat'); 
f_hco    = fopen('../HCO_Average_Scattering_Angles.dat'); 
f_hco2   = fopen('../HCO2_Average_Scattering_Angles.dat'); 
f_hear   = fopen('../HeAr_Average_Scattering_Angles.dat'); 
f_heh2   = fopen('../HeH2_Average_Scattering_Angles.dat'); 
f_hen2   = fopen('../HeN2_Average_Scattering_Angles.dat'); 
f_heco   = fopen('../HeCO_Average_Scattering_Angles.dat'); 
f_heco2  = fopen('../HeCO2_Average_Scattering_Angles.dat'); 

d_heo    = fscanf(f_heo,'%f %f',[2,inf]); 
d_hehe   = fscanf(f_hehe,'%f %f',[2,inf]); 
d_heh    = fscanf(f_heh,'%f %f',[2,inf]); 
d_hh     = fscanf(f_hh,'%f %f',[2,inf]); 
d_ho     = fscanf(f_ho,'%f %f',[2,inf]); 
d_har    = fscanf(f_har,'%f %f',[2,inf]); 
d_hh2    = fscanf(f_hh2,'%f %f',[2,inf]); 
d_hn2    = fscanf(f_hn2,'%f %f',[2,inf]); 
d_hco    = fscanf(f_hco,'%f %f',[2,inf]); 
d_hco2   = fscanf(f_hco2,'%f %f',[2,inf]); 
d_hear   = fscanf(f_hear,'%f %f',[2,inf]); 
d_heh2   = fscanf(f_heh2,'%f %f',[2,inf]); 
d_hen2   = fscanf(f_hen2,'%f %f',[2,inf]); 
d_heco   = fscanf(f_heco,'%f %f',[2,inf]); 
d_heco2  = fscanf(f_heco2,'%f %f',[2,inf]); 

mH 	 = 1; 
mHe  = 4; 
mO   = 16; 
mAr  = 40; 
mH2  = 2; 
mN2  = 28; 
mCO  = 28; 
mCO2 = 44; 

d_heo   = d_heo'; 
d_hehe  = d_hehe'; 
d_heh   = d_heh'; 
d_hh    = d_hh'; 
d_ho    = d_ho'; 
d_har    = d_har'; 
d_hh2    = d_hh2'; 
d_hn2    = d_hn2'; 
d_hco    = d_hco'; 
d_hco2    = d_hco2'; 
d_hear    = d_hear'; 
d_heh2    = d_heh2'; 
d_hen2    = d_hen2'; 
d_heco    = d_heco'; 
d_heco2    = d_heco2'; 

CC = pi/180; 

for i=1:length(d_heo)
	LE_heo(i) 	= d_heo(i,1)*(mHe+mO)/mO; 	
	LE_hehe(i) 	= d_hehe(i,1)*(mHe+mHe)/mHe; 	
	LE_heh(i) 	= d_heh(i,1)*(mHe+mH)/mH; 	
	LE_hhe(i) 	= d_heh(i,1)*(mHe+mH)/mHe; 	
	LE_hh(i) 		= d_hh(i,1)*(mH+mH)/mH; 	
	LE_ho(i) 		= d_ho(i,1)*(mH+mO)/mO; 	
	LE_har(i) 	= d_har(i,1)*(mH+mAr)/mAr; 	
	LE_hh2(i) 	= d_hh2(i,1)*(mH+mH2)/mH2; 	
	LE_hn2(i) 	= d_hn2(i,1)*(mH+mN2)/mN2; 	
	LE_hco(i) 	= d_hco(i,1)*(mH+mCO)/mCO; 	
	LE_hco2(i) 	= d_hco2(i,1)*(mH+mCO2)/mCO2; 	
	LE_hear(i) 	= d_hear(i,1)*(mHe+mAr)/mAr; 	
	LE_heh2(i) 	= d_heh2(i,1)*(mHe+mH2)/mH2; 	
	LE_hen2(i) 	= d_hen2(i,1)*(mHe+mN2)/mN2; 	
	LE_heco(i) 	= d_heco(i,1)*(mHe+mCO)/mCO; 	
	LE_heco2(i) = d_heco2(i,1)*(mHe+mCO2)/mCO2; 	

	t          	= d_heo(i,2)*CC; 
	p 				 	= mHe/mO; 
	LA_heo(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_hehe(i,2)*CC; 
	p 				 	= mHe/mHe; 
	LA_hehe(i) 	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_heh(i,2)*CC; 
	p 				 	= mHe/mH; 
	LA_heh(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_heh(i,2)*CC; 
	p 				 	= mH/mHe; 
	LA_hhe(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_hh(i,2)*CC; 
	p 				 	= mH/mH; 
	LA_hh(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_ho(i,2)*CC; 
	p 				 	= mH/mO; 
	LA_ho(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_har(i,2)*CC; 
	p 				 	= mH/mAr; 
	LA_har(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_hh2(i,2)*CC; 
	p 				 	= mH/mH2; 
	LA_hh2(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_hn2(i,2)*CC; 
	p 				 	= mH/mN2; 
	LA_hn2(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_hco(i,2)*CC; 
	p 				 	= mH/mCO; 
	LA_hco(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_hco2(i,2)*CC; 
	p 				 	= mH/mCO2; 
	LA_hco2(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_hear(i,2)*CC; 
	p 				 	= mHe/mAr; 
	LA_hear(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_heh2(i,2)*CC; 
	p 				 	= mHe/mH2; 
	LA_heh2(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_hen2(i,2)*CC; 
	p 				 	= mHe/mN2; 
	LA_hen2(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_heco(i,2)*CC; 
	p 				 	= mHe/mCO; 
	LA_heco(i)  	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
	t          	= d_heco2(i,2)*CC; 
	p 				 	= mHe/mCO2; 
	LA_heco2(i)	= acos((cos(t)+p)/sqrt(1+p^2+2*p*cos(t)))/CC; 
end

figure
h = semilogy( d_heo(:,1),d_heo(:,2),'k', ...
				  		d_hehe(:,1),d_hehe(:,2),'b', ...
				  		d_heh(:,1),d_heh(:,2),'g', ...
				  		d_hh(:,1),d_hh(:,2),'r', ... 
				  		d_ho(:,1),d_ho(:,2),'c' , ...
				  		d_har(:,1),d_har(:,2),'k--' , ...
				  		d_hh2(:,1),d_hh2(:,2),'b--' , ...
				  		d_hn2(:,1),d_hn2(:,2),'g--' , ...
				  		d_hco(:,1),d_hco(:,2),'r--' , ...
				  		d_hco2(:,1),d_hco2(:,2),'c--' , ...
				  		d_hear(:,1),d_hear(:,2),'k-.' , ...
				  		d_heh2(:,1),d_heh2(:,2),'b-.' , ...
				  		d_hen2(:,1),d_hen2(:,2),'g-.' , ...
				  		d_heco(:,1),d_heco(:,2),'r-.' , ...
				  		d_heco2(:,1),d_heco2(:,2),'c-.'); 
set(h,'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('CM Collision Energy [eV]')
ylabel('Average Scattering Angle [deg]')
legend('He+O', 'He+He', 'He+H', 'H+H','H+O','H+Ar',...
'H+H2','H+N2','H+CO','H+CO2','He+Ar','He+H2','He+N2',...
'He+CO','He+CO2','Location','Best')
print -depsc2 ./Plots/ALL_Average_Scattering_Angles.eps

figure
h = semilogy( LE_hh(:),  LA_hh(:),'ko', ...
				  		LE_hhe(:), LA_hhe(:),'bo', ... 
				  		LE_ho(:),  LA_ho(:),'go' , ...
				  		LE_har(:), LA_har(:),'ro' , ...
				  		LE_hh2(:), LA_hh2(:),'k+' , ...
				  		LE_hn2(:), LA_hn2(:),'b+' , ...
				  		LE_hco(:), LA_hco(:),'g+' , ...
				  		LE_hco2(:),LA_hco2(:),'r+','LineWidth',2.0); 
set(gca,'FontSize',16); 
xlabel('LAB CM Collision Energy [eV]')
ylabel('LAB Average Scattering Angle [deg]')
legend('H+H', 'H+He','H+O','H+Ar',...
'H+H2','H+N2','H+CO','H+CO2', 'Location','Best')
%axis([min(d_heh(:,1)) max(d_heh(:,1)) 5e-3 1e0])
print -depsc2 ./Plots/H_Average_Scattering_Angles_LAB.eps

figure
h = semilogy( LE_heh(:),  LA_heh(:),'ko', ...
				  		LE_hehe(:), LA_hehe(:),'bo', ... 
				  		LE_heo(:),  LA_heo(:),'go' , ...
				  		LE_hear(:), LA_hear(:),'ro' , ...
				  		LE_heh2(:), LA_heh2(:),'k+' , ...
				  		LE_hen2(:), LA_hen2(:),'b+' , ...
				  		LE_heco(:), LA_heco(:),'g+' , ...
				  		LE_heco2(:),LA_heco2(:),'r+','LineWidth',2.0); 
set(gca,'FontSize',16); 
xlabel('LAB CM Collision Energy [eV]')
ylabel('LAB Average Scattering Angle [deg]')
legend('He+H', 'He+He','He+O','He+Ar',...
'He+H2','He+N2','He+CO','He+CO2', 'Location','Best')
%axis([min(d_heh(:,1)) max(d_heh(:,1)) 5e-3 1e0])
print -depsc2 ./Plots/He_Average_Scattering_Angles_LAB.eps

figure
h = semilogy( d_hh(:,1),d_hh(:,2),'ko', ...
				  		d_heh(:,1),d_heh(:,2),'bo', ... 
				  		d_ho(:,1),d_ho(:,2),'go' , ...
				  		d_har(:,1),d_har(:,2),'ro' , ...
				  		d_hh2(:,1),d_hh2(:,2),'k+' , ...
				  		d_hn2(:,1),d_hn2(:,2),'b+' , ...
				  		d_hco(:,1),d_hco(:,2),'g+' , ...
				  		d_hco2(:,1),d_hco2(:,2),'r+','LineWidth',2.0); 
set(gca,'FontSize',16); 
xlabel('CM Collision Energy [eV]')
ylabel('CM Average Scattering Angle [deg]')
legend('H+H', 'H+He','H+O','H+Ar',...
'H+H2','H+N2','H+CO','H+CO2', 'Location','Best')
axis([min(d_heh(:,1)) max(d_heh(:,1)) 5e-3 1e0])
print -depsc2 ./Plots/H_Average_Scattering_Angles.eps

figure
h = semilogy( d_heh(:,1),d_heh(:,2),'ko', ...
				  		d_hehe(:,1),d_hehe(:,2),'bo', ...
				  		d_heo(:,1),d_heo(:,2),'go', ...
				  		d_hear(:,1),d_hear(:,2),'ro' , ...
				  		d_heh2(:,1),d_heh2(:,2),'k+' , ...
				  		d_hen2(:,1),d_hen2(:,2),'b+' , ...
				  		d_heco(:,1),d_heco(:,2),'g+' , ...
				  		d_heco2(:,1),d_heco2(:,2),'r+','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('CM Collision Energy [eV]')
ylabel('Average Scattering Angle [deg]')
legend('He+H', 'He+He', 'He+O', 'He+Ar',...
'He+H2','He+N2','He+CO','He+CO2','Location','Best')
axis([min(d_heh(:,1)) max(d_heh(:,1)) 5e-3 1e0])
print -depsc2 ./Plots/He_Average_Scattering_Angles.eps

mu1 = 4/5; 
mu2 = 4*4/8; 
mu3 = 4*16/20; 
mu4 = 4*40/44; 
mu5 = 4*2/6; 
mu6 = 4*28/32; 
mu7 = 4*28/32; 
mu8 = 4*44/48; 

figure
h = semilogy( sqrt(d_heh(:,1)/mu1/2),d_heh(:,2),'ko', ...
				  		sqrt(d_hehe(:,1)/mu2/2),d_hehe(:,2),'bo', ...
				  		sqrt(d_heo(:,1)/mu3/2),d_heo(:,2),'go', ...
				  		sqrt(d_hear(:,1)/mu4/2),d_hear(:,2),'ro' , ...
				  		sqrt(d_heh2(:,1)/mu5/2),d_heh2(:,2),'k+' , ...
				  		sqrt(d_hen2(:,1)/mu6/2),d_hen2(:,2),'b+' , ...
				  		sqrt(d_heco(:,1)/mu7/2),d_heco(:,2),'g+' , ...
				  		sqrt(d_heco2(:,1)/mu8/2),d_heco2(:,2),'r+','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('CM Collision Velocity [eV/amu]')
ylabel('Average Scattering Angle [deg]')
legend('He+H', 'He+He', 'He+O', 'He+Ar',...
'He+H2','He+N2','He+CO','He+CO2','Location','Best')
%axis([min(d_heh(:,1)) max(d_heh(:,1)) 5e-3 1e0])
print -depsc2 ./Plots/He_Average_Scattering_Angles.eps

mu1 = 4/5; 
mu2 = 4*16/20; 

figure
h = semilogy( d_heh(:,1)/mu1,d_heh(:,2),'ko', ...
				  		d_heo(:,1)/mu2,d_heo(:,2),'go', ...
							'LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('CM Collision Energy [eV]')
ylabel('Average Scattering Angle [deg]')
legend('He+H','He+O','Location','Best')
