
clear all
close all

f2 = fopen('../200km/planet_3d_Thermal_Height.dat'); 
d2 = fscanf(f2,'%f',[1,inf]); 
d2 = d2'; 

f3 = fopen('../300km/planet_3d_Thermal_Height.dat'); 
d3 = fscanf(f3,'%f',[1,inf]); 
d3 = d3'; 

f4 = fopen('../400km/planet_3d_Thermal_Height.dat'); 
d4 = fscanf(f4,'%f',[1,inf]); 
d4 = d4'; 

[y2,x2] = hist(d2,100); 
[y3,x3] = hist(d3,100); 
[y4,x4] = hist(d4,100); 

x2 = x2/1000; 
x3 = x3/1000; 
x4 = x4/1000; 

figure
h2 = plot(y2, x2,'k', ...
          y3, x3,'b', ...
          y4, x4,'r'); 
set(h2,'LineWidth',2.5)
xlabel('Frequency')
ylabel('Thermalization height [km]')
legend('200 km', '300 km', '400 km')


