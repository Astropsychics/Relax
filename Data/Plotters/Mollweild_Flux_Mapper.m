
clear all
close all

fid = fopen('../Flux_Map.dat'); 
dat = fscanf(fid,'%f %f %f %f',[4,inf]); 
dat = dat'; 

x = dat(:,1); 
y = dat(:,2); 
z = dat(:,3); 
E = dat(:,4); 
N = length(x); 
R(1:N) = 0; 
dif    = pi/1800; 

figure
for i=1:N
	R(i) 		 = sqrt(x(i)^2+y(i)^2+z(i)^2); 
	Theta(i) = acos(z(i)/R(i)); 
	Phi(i)   = atan(y(i)/x(i)); 
	P(i) = asin(y(i)/R(i)); 
	L(i) = asin(z(i)/cos(P(i))); 
	ddif = 1;  
	tn0  = P(i); 
	while (ddif > dif)
		tn1  = tn0 - ( 2*tn0+sin(2*tn0)-pi*sin(P(i)))/(2+2*cos(2*tn0)); 
		ddif = abs(tn1-tn0); 
		tn0  = tn1;  	
	end
	T(i) = tn0; 
	X(i) = (2*sqrt(2)*L(i)/pi)*cos(T(i)); 
	Y(i) = sqrt(2)*sin(T(i)); 

	if ( E(i) <= 500 )
		plot(X(i),Y(i),'b*','LineWidth',2.5) 
	elseif ( E(i) <= 1000 )
		plot(X(i),Y(i),'g*','LineWidth',2.5) 
	elseif ( E(i) <= 1500 )
		plot(X(i),Y(i),'c*','LineWidth',2.5) 
	else 
		plot(X(i),Y(i),'r*','LineWidth',2.5) 
	end
	hold on
	
end

set(gca,'FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
print -depsc2 ./Plots/Flux.eps

figure
for i=1:N
  if ( E(i) <= 500 )
    plot3(x(i),y(i),z(i),'b*','LineWidth',2.5)
  elseif ( E(i) <= 1000 )
    plot3(x(i),y(i),z(i),'g*','LineWidth',2.5)
  elseif ( E(i) <= 1500 )
    plot3(x(i),y(i),z(i),'c*','LineWidth',2.5)
  else
    plot3(x(i),y(i),z(i),'r*','LineWidth',2.5)
  end
  hold on
 
end

figure
for i=1:N
  if ( E(i) <= 500 )
    plot(Theta(i),Phi(i),'b*','LineWidth',2.5)
  elseif ( E(i) <= 1000 )
    plot(Theta(i),Phi(i),'g*','LineWidth',2.5)
  elseif ( E(i) <= 1500 )
    plot(Theta(i),Phi(i),'c*','LineWidth',2.5)
  else
    plot(Theta(i),Phi(i),'r*','LineWidth',2.5)
  end
  hold on
end



