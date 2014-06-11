
clear all
close all

fid = fopen('../planet_3d_Collision_Number.dat'); 
dat = fscanf(fid,'%f',[1,inf]); 
dat = dat'; 

figure
hist(dat, 50)


