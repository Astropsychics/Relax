
clear all
close all

f = fopen('../Escape_Energy_Distribution.dat'); 
d = fscanf(f,'%f',[1,inf]); 
d = d'; 

NH = 20; 

[y,x] = hist(d,NH); 

dx = x(2)-x(1); 
N  = length(d); 
y  = y/(sum(y)*dx); 

figure
plot(x,y,'k',x,y,'go','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Escape Energy [eV]')
ylabel('Probability [1/eV]')
print -depsc2 ./Plots/Escape_Energy_PDF.eps

for i=1:NH
	tot = 0; 
	for j=1:i
		tot = tot + y(j)*dx; 
	end
	cum_y(i) = tot; 
end

figure
plot(cum_y, x,'k',cum_y,x,'ro','LineWidth',2.5)
set(gca,'FontSize',16)
ylabel('Escape Energy [eV]')
xlabel('Cumulative Probability')
axis([0 1 min(x) max(x)])
print -depsc2 ./Plots/Escape_Energy_CDF.eps

	


