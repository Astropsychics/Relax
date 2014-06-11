
clear all
close all

xN  = '../planet_X_Ncoll.dat';
xE  = '../planet_X_energy.dat';
xdE = '../planet_X_energy_loss.dat';
xH  = '../planet_X_height.dat';
xP  = '../planet_X_phi.dat';
xA  = '../planet_X_theta.dat';
xT  = '../planet_X_time.dat';
xU  = '../planet_X_unit_vel.dat'; 

f1 = '../planet_all_Ncoll_vs_time.dat'; 
x1 = xT; 
y1 = xN; 
N1 = '<N> vs T'; 
f2 = '../planet_all_energy_vs_time.dat';
x2 = xT; 
y2 = xE; 
N2 = '<E> vs T'; 
f3 = '../planet_all_height_vs_Ncoll.dat';
x3 = xH; 
y3 = xN; 
N3 = '<N> vs Z'; 
f4 = '../planet_all_height_vs_energy.dat';
x4 = xH; 
y4 = xE; 
N4 = '<E> vs Z'; 
f5 = '../planet_all_height_vs_energy_loss.dat';
x5 = xH; 
y5 = xdE; 
N5 = '<dE> vs Z'; 
f6 = '../planet_all_height_vs_time.dat';
x6 = xT; 
y6 = xH; 
N6 = '<H> vs T'; 
f7 = '../planet_all_height_vs_vert_vel.dat';
x7 = xH; 
y7 = xU; 
N7 = '<Uz> vs Z'; 
f8 = '../planet_all_vert_vel_vs_time.dat'; 
x8 = xT; 
y8 = xU; 
N8 = '<Uz> vs T'; 

Click_ALL_Gen(f1,x1,y1,N1)
Click_ALL_Gen(f2,x2,y2,N2)
Click_ALL_Gen(f3,x3,y3,N3)
Click_ALL_Gen(f4,x4,y4,N4)
Click_ALL_Gen(f5,x5,y5,N5)
Click_ALL_Gen(f6,x6,y6,N6)
Click_ALL_Gen(f7,x7,y7,N7)
Click_ALL_Gen(f8,x8,y8,N8)



