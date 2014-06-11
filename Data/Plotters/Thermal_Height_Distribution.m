
clear all
close all

f = fopen('../Thermalization_Height_Distribution.dat'); 
d = fscanf(f,'%f',[1,inf]); 
d = d'; 

NH = 20; 

[y,x] = hist(d,NH); 

N = length(d); 
y = 100*y/N; 
x = x/1000; 

figure
plot(x,y,'k',x,y,'go','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Thermalization Height [km]')
ylabel('Ensemble Percentage')
print -depsc2 ./Plots/Thermal_Height_PDF.eps

for i=1:NH
	tot = 0; 
	for j=1:i
		tot = tot + y(j); 
	end
	cum_y(i) = tot/100; 
end

figure
plot(cum_y, x,'k',cum_y,x,'ro','LineWidth',2.5)
set(gca,'FontSize',16)
ylabel('Thermalization Height [km]')
xlabel('Cumulative Probability')
axis([0 1 min(x) max(x)])
print -depsc2 ./Plots/Thermal_Height_CDF.eps

	


