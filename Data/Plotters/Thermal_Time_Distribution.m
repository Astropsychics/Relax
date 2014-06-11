
clear all
close all

f = fopen('../Thermalization_Time_Distribution.dat'); 
d = fscanf(f,'%f',[1,inf]); 
d = d'; 

NH = 20; 

[y,x] = hist(d,NH); 

N = length(d); 
y = 100*y/N; 

figure
plot(x,y,'k',x,y,'go','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Thermalization Time [sec]')
ylabel('Ensemble Percentage')
print -depsc2 ./Plots/Thermal_Time_PDF.eps

for i=1:NH
	tot = 0; 
	for j=1:i
		tot = tot + y(j); 
	end
	cum_y(i) = tot/100; 
end

X(1) = 0; 
X(2:length(x)+1) = x; 
C(1) = 0; 
C(2:length(x)+1) = cum_y; 

figure
%plot(cum_y, x,'k',cum_y,x,'ro','LineWidth',2.5)
plot(C,X,'k',C,X,'ro','LineWidth',2.5)
set(gca,'FontSize',16)
ylabel('Thermalization Time [sec]')
xlabel('Cumulative Probability')
axis([0 1 min(X) max(X)])
print -depsc2 ./Plots/Thermal_Time_CDF.eps

	


