close all;clear all;clc

Data_in     =   './forepaw_36_perts_loc_01';
Jinv_in     =   './jinv_file';
Sol_out     =   './Forepaw_';

load(Jinv_in,'Jinv')
load(Data_in,'bnd_v','bnd_v_c0','bnd_v_c')
tic
for kk=1%:size(bnd_v_c,2)
    
    diff_data0  =	100*( bnd_v_c0(:,kk) - bnd_v)./bnd_v;
    diff_data   =   100*( bnd_v_c(:,kk) - bnd_v)./bnd_v;
    
    sol0        =   -Jinv*diff_data0; %noise free %5 mins
    sol         =   -Jinv*diff_data;% with noise
    
    save([Sol_out 'pos_' num2str(kk)],'sol*')
    
    clear sol*
    
end
toc