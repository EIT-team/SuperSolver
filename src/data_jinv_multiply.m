close all;clear all;clc

Data_in     =   'D:\Dropbox\PhD\Research\Neural Imaging\03_Fast_Neural_Trajectories\04_Simulation\03_Images\01_BVs\forepaw_36_perts_loc_01';
Jinv_in     =   'D:\Neural Imaging\03_Fast_Neural_Trajectories\04_Simulation\01_Inverse_Crime\02_Jinvs_and_Sols\Adap_Mesh50';
Sol_out     =   'D:\Neural Imaging\03_Fast_Neural_Trajectories\04_Simulation\03_Images\01_Solutions\Forepaw_';

load(Jinv_in,'Jinv')
load(Data_in,'bnd_v','bnd_v_c0','bnd_v_c')

for kk=5:size(bnd_v_c,2)
    
    diff_data0  =	100*( bnd_v_c0(:,kk) - bnd_v)./bnd_v;
    diff_data   =   100*( bnd_v_c(:,kk) - bnd_v)./bnd_v;
    
    sol0        =   -Jinv*diff_data0; %noise free %5 mins
    sol         =   -Jinv*diff_data;% with noise
    
    save([Sol_out 'pos_' num2str(kk)],'sol*')
    
    clear sol*
    
end