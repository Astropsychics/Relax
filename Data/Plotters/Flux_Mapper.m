
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

for i=1:N
	R(i) 		 = sqrt(x(i)^2+y(i)^2+z(i)^2); 
	Theta(i) = acos(z(i)/R(i)); 
	Phi(i)   = atan(y(i)/x(i)); 
end


figure
for i=1:N
  if ( E(i) <= 500 )
    plot(Phi(i),Theta(i),'b*','LineWidth',2.5)
  elseif ( E(i) <= 1000 )
    plot(Phi(i),Theta(i),'g*','LineWidth',2.5)
  elseif ( E(i) <= 1500 )
    plot(Phi(i),Theta(i),'c*','LineWidth',2.5)
  else
    plot(Phi(i),Theta(i),'r*','LineWidth',2.5)
  end
  hold on
end
set(gca,'FontSize',16)
xlabel('Phi [rad]')
ylabel('Theta [rad]')
axis([-pi/2 pi/2 0 pi])


N_Grid = 20; 
M_Grid = N_Grid - 1; 
t_in   = -pi/2; 
t_fn   = pi/2; 
p_in   = -pi/2; 
p_fn   = pi/2; 
dtheta = (t_fn-t_in)/M_Grid; 
dphi   = (p_fn-p_in)/M_Grid; 

for i=1:N_Grid
	G_Theta(i) = t_in + (i-1)*dtheta; 
	G_Phi(i)   = p_in + (i-1)*dphi; 
end

G_E(1:M_Grid,1:M_Grid) = 0; 
G_C(1:M_Grid,1:M_Grid) = 0; 

for i=1:N
	Pnow = Phi(i); 
	Tnow = Theta(i); 
	for j=1:M_Grid
		P1 = G_Phi(j); 
		P2 = G_Phi(j+1); 
		for k=1:M_Grid
			T1 = G_Theta(k); 
			T2 = G_Theta(k+1); 
			if ( ((P1 <= Pnow) && (Pnow <= P2)) && ((T1 <= Tnow) && (Tnow <= T2)) )
				G_E(j,k) = G_E(j,k) + E(i); 	
				G_C(j,k) = G_C(j,k) + 1; 	
			end		
		end
	end
end

for i=1:M_Grid
	for j=1:M_Grid
		if (G_C(i,j) ~= 0)
			FLUX(i,j) = G_E(i,j)/G_C(i,j); 
		else
			Flux(i,j) = 0; 
		end
	end
end

figure
imagesc(G_Phi,G_Theta,FLUX')
set(gca,'FontSize',16,'YDir','Normal')
xlabel('Phi [rad]')
ylabel('Theta [rad]')
colorbar





