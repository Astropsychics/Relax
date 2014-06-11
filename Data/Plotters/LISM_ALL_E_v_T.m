
clear all
close all

xN  = '../lism_X_Ncoll.dat';
xE  = '../lism_X_energy.dat';
xdE = '../lism_X_energy_loss.dat';
xH  = '../lism_X_height.dat';
xP  = '../lism_X_phi.dat';
xA  = '../lism_X_theta.dat';
xT  = '../lism_X_time.dat';
xU  = '../lism_X_unit_vel.dat'; 

W1 = 0; 
W2 = 1; 
W3 = 0; 
W4 = 0; 
W5 = 0; 
W6 = 0; 
W7 = 0; 
W8 = 0; 

if (W1 == 1) 
f1 = '../lism_all_Ncoll_vs_time.dat'; 
x1 = xT; 
y1 = xN; 
N1 = '<N> vs T'; 
Click_ALL_Gen(f1,x1,y1,N1)
end

if (W2 == 1) 
f2 = '../lism_all_energy_vs_time.dat';
x2 = xT; 
y2 = xE; 
N2 = '<E> vs T'; 
Click_ALL_Gen(f2,x2,y2,N2)
end

if (W3 == 1) 
f3 = '../lism_all_height_vs_Ncoll.dat';
x3 = xH; 
y3 = xN; 
N3 = '<N> vs Z'; 
Click_ALL_Gen(f3,x3,y3,N3)
end

if (W4 == 1) 
f4 = '../lism_all_height_vs_energy.dat';
x4 = xH; 
y4 = xE; 
N4 = '<E> vs Z'; 
Click_ALL_Gen(f4,x4,y4,N4)
end

if (W5 == 1) 
f5 = '../lism_all_height_vs_energy_loss.dat';
x5 = xH; 
y5 = xdE; 
N5 = '<dE> vs Z'; 
Click_ALL_Gen(f5,x5,y5,N5)
end

if (W6 == 1) 
f6 = '../lism_all_height_vs_time.dat';
x6 = xT; 
y6 = xH; 
N6 = '<H> vs T'; 
Click_ALL_Gen(f6,x6,y6,N6)
end

if (W7 == 1) 
f7 = '../lism_all_height_vs_vert_vel.dat';
x7 = xH; 
y7 = xU; 
N7 = '<Uz> vs Z'; 
Click_ALL_Gen(f7,x7,y7,N7)
end

if (W8 == 1) 
f8 = '../lism_all_vert_vel_vs_time.dat'; 
x8 = xT; 
y8 = xU; 
N8 = '<Uz> vs T'; 
Click_ALL_Gen(f8,x8,y8,N8)
end




