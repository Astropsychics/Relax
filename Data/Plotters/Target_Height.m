
clear all
close all

fid = fopen('../planet_3d_Target_Height.dat'); 
dat = fscanf(fid,'%f %f',[2,inf]); 
dat = dat'; 

iH   = 1; 
iHe  = 1; 
iO   = 1; 
iCO2 = 1; 

for i=1:length(dat(:,1))

	if ( dat(i,2) == 1 )
		H(iH) = dat(i,1); 
		iH = iH + 1; 

	elseif ( dat(i,2) == 4 )
		He(iHe) = dat(i,1); 
		iHe = iHe + 1; 

	elseif ( dat(i,2) == 16 ) 
		O(iO) = dat(i,1); 
		iO = iO + 1; 

	elseif ( dat(i,2) == 44 )
		CO2(iCO2) = dat(i,1); 
		iCO2 = iCO2 + 1; 

	end

end


