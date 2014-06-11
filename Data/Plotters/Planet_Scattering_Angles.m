
clear all
close all

f_co2 = fopen('../planet_CO2_Scatt_Angle.dat'); 
f_o   = fopen('../planet_O_Scatt_Angle.dat'); 
f_he  = fopen('../planet_He_Scatt_Angle.dat'); 
f_h   = fopen('../planet_H_Scatt_Angle.dat'); 

d_co2 = fscanf(f_co2,'%f %f',[2,inf]); 
d_o   = fscanf(f_o,'%f %f',[2,inf]); 
d_he  = fscanf(f_he,'%f %f',[2,inf]); 
d_h   = fscanf(f_h,'%f %f',[2,inf]); 

d_co2 = d_co2'; 
d_o   = d_o'; 
d_he  = d_he'; 
d_h   = d_h'; 

co2_max = max(d_co2(:,1)); 
o_max   = max(d_o(:,1)); 
he_max  = max(d_he(:,1)); 
h_max   = max(d_h(:,1)); 

En_N   = 51; 

dE_co2 = co2_max/En_N; 
dE_o   = o_max/En_N; 
dE_he  = he_max/En_N; 
dE_h   = h_max/En_N; 

for i=1:(En_N-1)
	E1_co2 = i*dE_co2; 
	E2_co2 = (i+1)*dE_co2; 
	E1_o   = i*dE_o; 
	E2_o   = (i+1)*dE_o; 
	E1_he  = i*dE_he; 
	E2_he  = (i+1)*dE_he; 
	E1_h   = i*dE_h; 
	E2_h   = (i+1)*dE_h; 

	E_co2(i) = E1_co2 + (E2_co2 - E1_co2)/2; 	
	E_o(i)   = E1_o   + (E2_o   - E1_o)  /2; 	
	E_he(i)  = E1_he  + (E2_he  - E1_he) /2; 	
	E_h(i)   = E1_h   + (E2_h   - E1_h)  /2; 	

%% CO2
	c = 0; 
	t = 0; 
	clear X_co2; 
	for j=1:length(d_co2(:,1))
		if ( (E1_co2 <= d_co2(j,1)) && (d_co2(j,1) <= E2_co2) )
			t = t + d_co2(j,2); 
			c = c + 1;
			X_co2(c) = d_co2(j,2);  
		end		
	end
	if (c > 0)
		Ang_co2(i) = t/c; 
	else
		Ang_co2(i) = 0; 
	end
	Std_co2(i) = std(X_co2); 

%% O
	c = 0; 
	t = 0; 
	clear X_o; 
	for j=1:length(d_o(:,1))
		if ( (E1_o <= d_o(j,1)) && (d_o(j,1) <= E2_o) )
			t = t + d_o(j,2); 
			c = c + 1;
			X_o(c) = d_o(j,2);  
		end		
	end
	if (c > 0)
		Ang_o(i) = t/c; 
	else
		Ang_o(i) = 0; 
	end
	Std_o(i) = std(X_o); 

%% He
	c = 0; 
	t = 0; 
	for j=1:length(d_he(:,1))
		if ( (E1_he <= d_he(j,1)) && (d_he(j,1) <= E2_he) )
			t = t + d_he(j,2); 
			c = c + 1; 
		end		
	end
	if (c > 0)
		Ang_he(i) = t/c; 
	else
		Ang_he(i) = 0; 
	end

%% H
	c = 0; 
	t = 0; 
	for j=1:length(d_h(:,1))
		if ( (E1_h <= d_h(j,1)) && (d_h(j,1) <= E2_h) )
			t = t + d_h(j,2); 
			c = c + 1; 
		end		
	end
	if (c > 0)
		Ang_h(i) = t/c; 
	else
		Ang_h(i) = 0; 
	end

end

figure
h = plot(E_co2,Ang_co2,'k', ...
         E_o,  Ang_o,  'b', ...
         E_he, Ang_he, 'r', ...
         E_h,  Ang_h,  'g'); 
set(h,'LineWidth',2.5); 
xlabel('Lab Collision Energy [eV]')
ylabel('Average Scattering Angle [deg]')
legend('CO2', 'O', 'He', 'H','Location','NorthEast')
print -depsc2 Average_Scattering_Angles.eps

figure
errorbar(E_co2, Ang_co2, Std_co2,'k')
hold on
errorbar(E_o, Ang_o, Std_o,'b')
xlabel('Lab Collision Energy [eV]')
ylabel('Average Scattering Angle [deg]')
legend('CO2','O')
print -depsc2 Average_Scattering_Angles_Errorbars.eps


