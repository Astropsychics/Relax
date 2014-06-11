
clear all
close all

f = fopen('../Initial_Energy.dat'); 
d = fscanf(f,'%f',[1,inf]); 
d = d'; 

NH = 100; 

[y,x] = hist(d,NH); 

N = length(d); 
y = 100*y/N; 

figure
plot(x,y,'go',x,y,'k','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Energy [eV]')
ylabel('Ensemble Percentage')
print -depsc2 ./Plots/Initital_Energy_PDF.eps

for i=1:NH
	tot = 0; 
	for j=1:i
		tot = tot + y(j); 
	end
	cum_y(i) = tot/100; 
end

figure
plot(x,cum_y,'ro',x,cum_y,'k','LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Energy [eV]')
ylabel('Cumulative Probability')
axis([min(x) max(x) 0 1])
print -depsc2 ./Plots/Initital_Energy_CDF.eps

	


