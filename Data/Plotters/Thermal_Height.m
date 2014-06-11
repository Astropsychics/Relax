
clear all
close all

PLOT_ON = 1; 

f = fopen('../planet_3d_Thermal_Height.dat'); 
d = fscanf(f,'%f',[1,inf]); 
d = d'; 


[freq,height] = hist(d,100); 

height = height/1000; 

figure
h = plot(freq, height,'k'); 
set(h,'LineWidth',2.5)
xlabel('Frequency')
ylabel('Thermalization height [km]')
if (PLOT_ON == 1) 
	print -depsc2 ./Plots/Thermalization_Heights.eps
end

