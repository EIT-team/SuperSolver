%clear all;close all;%clc
clear Mesh;
disp('Setup started')
tic
%==========================================================================%
% SWITCH SETUP SECTION %
LINUX       =       0;  % 0 - windows, 1 - linux
MODEL       =       3;  % 1 - cube, 2 - human head, 3 - rat brain
STIM        =       1;  % 1 - forepaw, 2 - whisker, 3 - hindpaw, 4 - visual
POS_29_31   =       1;  % 1 - 29 array; 0 - 31 array
SAVE        =       0;
ADD_NOISE   =       0;
Noise_Level =       .1; % noise level in uV
EXP_DATA    =       0;  % 0 - simulated data with 962 protocol,
% 1 - also solves for experimental data in addition to the default simulated data but with the experimental protocol
SIM         =       0;  % 0 - no simulation, 1 - simulations
BATCH_NO    =       1;
svd_list    =       [50];
VOLUME      =       0;  % 0 - finds a single point/tetrahderon, 1 - finds a blob of tetrahedra
rec         =       36;
FUCK_YOU    =       1;

SOLVE       =       1;  % 0 - stops after boundary voltage generation, 1 - runs through to inversion

%==========================================================================%
% PATH SETUP SECTION %
% Input_Files_Path='D:\Dropbox\PhD\Research\Neural Imaging\03_Fast_Neural_Trajectories\04_Simulation\00_Code\01_Input_Variables\';
% Output_Data_Path='D:\Dropbox\PhD\Research\Neural Imaging\03_Fast_Neural_Trajectories\04_Simulation\00_Code\02_Output_Variables\BVs\';
% Output_Sol_Path ='D:\Neural Imaging\03_Fast_Neural_Trajectories\04_Simulation\03_Images\01_Solutions\Rec36_Prt\Rec36_tet1';
tag             =   'forepaw_36_perts_loc_01';

if LINUX
    %     cd('./SimSolver/Code/')
    Input_Files_Path   =    linux_path(Input_Files_Path);
    Output_Data_Path   =    linux_path(Output_Data_Path);
    Output_Sol_Path    =    linux_path(Output_Sol_Path);
end

%==========================================================================%
% DEFINE MESH PARAMETERS %
global Mesh
Mesh.fn                             =   'Ilan_Rat_Brain_0p2to0p8mm_3V_1072558el';
Mesh.vtx                            =   []; % =vtx;
Mesh.units                          =   'm'; % meters
Mesh.tri                            =   []; %=tri;
Mesh.srf                            =   []; % =srf;
Mesh.cnts                           =   []; % =cnts; % elemetnts centres
Mesh.maps                           =   'element'; % 'nodal' 'regular'
Mesh.hgrid                          =   []; % =hgrid;
Mesh.rgrid.xi                       =   []; % =xi;
Mesh.rgrid.yi                       =   []; % =yi;
Mesh.rgrid.zi                       =   []; % =zi;
Mesh.sgrad                          =   []; % =D; % shape gradient functions
Mesh.volumes                        =   []; % elements volumes
Mesh.reg                            =   []; % =Reg; Gaussian regularisation matrix

% bla=load('.\mesh_prot_elec\Mesh_refined_electrode_30.mat');
% %bla=load('.\mesh_prot_elec\Mesh_refined_all_electrodes_coarse.mat');
% %bla=load('.\mesh_prot_elec\Mesh_refined_all_electrodes_slightly_finer_with_bigger_radius.mat');
% %bla=load('.\mesh_prot_elec\Mesh_refined_all_electrodes_fine_with_smaller_radius.mat');
% %bla = load('.\mesh_prot_elec\Mesh_huge_coarsely_refined_electrode_30.mat');
% %bla = load('.\mesh_prot_elec\Mesh_refined_all_electrodes_coarse_with_big_radius.mat');
% %bla = load('.\mesh_prot_elec\Mesh_refined_all_electrodes_even_coarser_with_large_radius.mat');
% Mesh.vtx = bla.Nodes(:,1:3)/1000;
% Mesh.tri = double(bla.Tetra(1:length(bla.Tetra)-length(find(bla.Tetra==0))/5,1:4)); %cut away zeros at the end of Tetra
% [Mesh.vtx, Mesh.tri] = removeisolatednode(Mesh.vtx, Mesh.tri);
% 
% mat_ref = double(bla.Tetra(1:length(bla.Tetra)-length(find(bla.Tetra==0))/5,5));
% mat_ref(mat_ref == 1) = 0.15;   % white matter
% mat_ref(mat_ref == 2) = 0.3;    % grey matter
% mat_ref(mat_ref == 3) = 1.79;   % CSF
% mat_ref(mat_ref == 4) = 0.44;   % dura mater
% mat_ref(mat_ref == 5) = 0.018;  % skull
% mat_ref(mat_ref == 6) = 0.0001; % air
% mat_ref(mat_ref == 7) = 0.44;   % white matter

bla = load('.\mesh_prot_elec\TA052.mat');
Mesh.vtx = bla.vtx/1000;
Mesh.tri = bla.tri;
mat_ref = bla.mat_ref;
pos = bla.pos/1000;

%mat_ref = 0.2*ones(size(mat_ref));

% bla = load('.\mesh_prot_elec\TA052.mat');
% Mesh.vtx = bla.vtx/1000.;
% Mesh.tri = bla.tri;
% mat_ref = bla.mat_ref;
% pos = bla.pos/1000.;

clear bla

%==========================================================================%
% DEFINE FEM PARAMETERS
global Fem
Fem.current                         =   133e-6; % 50e-6;% [50uA in rat and tank expt, 400uA is some simulations]
% Fem.elec_diam                       =   .6e-3;% electrode diamter in meters
Fem.elec_diam                       =   7e-3;% electrode diameter in meters
Fem.materials_fn                    =   [];
Fem.elec_method                     =   's'; % electrodes assigning method 'n', 's', 'm', 'v' closest node, closest elemt centre, nodes within Elec_diam/2 from remote-most surface coordinate, nodes within Elec_diam/2
Fem.df                              =   []; % df;
Fem.elec                            =   []; % elec;
Fem.gnd                             =   []; % gnd_ind;
Fem.I                               =   []; % I - sources matrix
Fem.zc                              =   1e3 ; % contact impedance (could be a vector at the length of the number of electrodes);
Fem.E                               =   []; % E; FEM system matrix
Fem.isotropy                        =   true; % define whether a model is Isotropic or anisotropic

% ==========================================================================%
% DEFINE THE PROTOCOL %
% if EXP_DATA
%     Fem.prt_fn                      =   'protocol_rem'; % 461 protocol
% else
%     Fem.prt_fn                      =   'Brick_prt'; % 962 protocol
% end

% load([Input_Files_Path Fem.prt_fn])

% load('D:\Dropbox\PhD\Research\Neural Imaging\03_Fast_Neural_Trajectories\04_Simulation\00_Code\01_Input_Variables\All_Rats_Imaging_Data.mat','Out_All')
% prt                                 =   Out_All{rec}.prt;
% clear Out_All

%Fem.prt                             =   prt;
%Fem.prt = [1 30 2 7];%; 1 2 3 4;1 3 2 4;1 4 2 3;2 3 1 4;2 4 1 3;3 4 1 2];% load('mesh_prot_elec/eeg31b.prt');%
%Fem.prt = load('mesh_prot_elec/eeg31b.prt');

% Fem.prt = [1	30	2	7   
% 1	30	7	17
% 1	30	17	27
% 1	30	29	23
% 1	30	23	13
% 1	30	13	3
% 1	30	18	31
% 1	30	31	28
% 1	30	28	12
% 1	30	4	9]; %POLAR INJECTION
% Fem.prt = [26	30	2	7   
% 26	30	7	17
% 26	30	17	27
% 26	30	29	23
% 26	30	23	13
% 26	30	13	3
% 26	30	18	31
% 26	30	31	28
% 26	30	28	12
% 26	30	4	9]; %ADJACENT INJECTION
% Fem.prt = [31	1	2	7   
% 31	1	7	17
% 31	1	17	27
% 31	1	29	23
% 31	1	23	13
% 31	1	13	3
% 31	1	18	30
% 31	1	30	28
% 31	1	28	12
% 31	1	4	9]; %Check electrode 31
%Fem.prt = [Fem.prt(:,3:4) Fem.prt(:,1:2)];
Fem.prt = [1, 30, 2, 7
    1, 2, 3, 4
    1, 3, 4, 2
    1, 4, 2, 3
    2, 3, 1, 4
    2, 4, 1, 3
    30, 1, 7, 2]; % TA052 test protocol

%==========================================================================%

% READ ELECTRODES POSITION FILE %
if MODEL==3
    switch STIM
        case 1 % forepaw
            Fem.pos_fn              =   'forepaw_pos';
        case 2 % whisker
            Fem.pos_fn              =   'whisker_pos';
        case 3 % hindpaw
            Fem.pos_fn              =   'hinpaw_pos';
        case 4 % visual
            Fem.pos_fn              =   'visual_pos';
    end
end


%% load('.\mesh_prot_elec\correct_pos.mat')
Fem.pos                             =   pos;
Mesh.pos                            =   pos;
Fem.gnd_pos                         =   [0, -0.104467562567864, -0.0480000000021847];%[0.1337000, 0.0125200, 0.0550700];%[0, -0.104467562567864, -0.0480000000021847];% %arbitrarily chosen as the CLA marker
Fem.elec_diam                       =   repmat(Fem.elec_diam,length(Fem.pos),1); % electrode diameter in meters as a vector
Fem.elec_diam(30)                   =   7.0e-3;
%Fem.elec_diam(30)                   =   15.00e-3 - 2.86e-3 - 1.32e-3 - 0.67e-3 USW; %adjacent injection. step size is 10 to increase convergence speed.

Sigma1_change_factor                =   1.01; % 1% simulated perturbation change
Sigma2_change_factor                =   1.00;

%==========================================================================%
% CHANGE POSITION OF ELECTRODE 30
%size_of_change = -70e-4 + 11.32e-4 + 12.26e-4 + 19.27e-4 + 9.72e-4 + 6.59e-4 + 8.07e-4 + 5.52e-4 + 8.45e-4 + 4.47e-4 + 4.16e-4 + 4.33e-4 + 5.81e-4 + 2.44e-4...
%                + 2.44e-4 + 3.22e-4 + 3.22e-4 +  3.11e-4 + 2.78e-4; % X these are with filtered Jacobian jac(<1e-5) = inf
%size_of_change = size_of_change + 2.78e-3 - 0.47e-3; % X from here on I used a step width of 1e-2 for faster convergence. ( i.e. (1./test(:,30))'*diff * -1e-3 )
%size_of_change = -70e-4 + 1.317e-3 + 1.631e-3 + 0.676e-3 + 1.236e-3 + 1.336e-3 + 0.521e-3 + 0.811e-3 + 0.409e-3; % Y stepsize is 1e-2 => inv(elec_jac(:,30)'*elec_jac(:,30)) * elec_jac(:,30)' * diff * -1e-2
%size_of_change = size_of_change + 7.722e-3 - 1.675e-3; % Y from here on I used a step width of 1e-1 for faster convergence.
%size_of_change = -70e-4 + 5.75e-4 + 45.83e-4 - 6.09e-4 - 5*4.15e-4 + 68.02e-4 - 3*3.07e-4 - 3*1.69e-4 - 5*4.98e-4; %this is for electrode 1 (very coarse). stepsize is 1e-2
%size_of_change = -70e-4 + 181.0e-4 - 56.09e-4 + 23.36e-4 - 15.02e-4 USW; %adjacent injection, step size is 1e-1
%size_of_change = 40e-4 - 17e-4 - 17e-4 - 8.5e-4 - 9.5e-4 - 5.3e-4 - 5.1e-4 - 5.1e-4 + -5e-4 - 5e-4; % electrode 1 which is very coarsely refined. Step size is 1e-2. Towards the end it goes in wrong direction at some point around -27e-4
%size_of_change = 40e-4 - 15.5e-4 - 20.7e-4 - 7.3e-4 - 8.2e-4 - 4.9e-4 - 4.0e-4 - 4.0e-4 - 5*1.2e-4 - 5*0.9e-4; %this is for electrode one using the full protocol.
size_of_change = 0e-4;
%size_of_change = -70e-4 + 202.8e-4 - 107.4e-4 + 47.8e-4 - 2e-4;
electrode_nr = 30;

[m,n]=size(coords);
coordinates = zeros(m*n/3,3);
for i = 1:n/3
    coordinates(1+(i-1)*m:i*m,:) = coords(:,1+(i-1)*3:i*3);
end

dist=(sum((coordinates-repmat(Fem.pos(electrode_nr,:),size(coordinates,1),1)).^2,2)).^0.5;
[~,idx]=min(dist);

areaORposition = 1;
if areaORposition == 1 % jacobian with respect to x-coordinate
    if dist(idx+m)<dist(idx-m) % electrode is closer to x+1 than x-1
        vector = coordinates(idx+m,:)-coordinates(idx,:);
        vector = vector/norm(vector)*size_of_change;
    else
        vector = coordinates(idx,:)-coordinates(idx-m,:);
        vector = vector/norm(vector)*size_of_change;
    end
else % jacobian with respect to y-coordinate
    if dist(idx+1)<dist(idx-1) % electrode is closer to y+1 than y-1
        vector = coordinates(idx+1,:)-coordinates(idx,:);
        vector = vector/norm(vector)*size_of_change;
    else
        vector = coordinates(idx,:)-coordinates(idx-1,:);
        vector = vector/norm(vector)*size_of_change;
    end
end

Fem.pos(electrode_nr,:) = Fem.pos(electrode_nr,:) + vector;


%==========================================================================%
% DEFINE FORWARD SOLVER PARAMETERS
global Fwd
Fwd.fn                              =   []; % file containing all the forward model parameters
Fwd.method                          =   'bicgstab'; % set forward solver method (either 'pcg', 'bicgstab', 'minres' or 'gmres' for pcg, stab BiCG, MinRes, GMRes)
Fwd.tol                             =   1e-12; % fwd solution relative residual tolerance (default 1e-12)
Fwd.maxit                           =   10000;%100; % forward solution maximum number of iterations (default 100)
Fwd.restarts                        =   20; % define restatrs for GMRes solver only (default 20)
Fwd.precond                         =   'ilu'; % set preconditioner method (either 'chol', 'amg', 'ilu', 'camg' for incomple Cholesky, AMG, incomplete LU, complex AMG)
Fwd.drop_tol                        =   1e-5; % drop tolerance fot incomplete preconditioners (default 1e-5)
Fwd.current_field                   =   []; % =V; % nodal current field potential distribution
Fwd.measurement_field               =   []; % =v_f; % nodal measurement field potential distribution

%==========================================================================%
% DEFINE SENSITIVITY PARAMETERS
global Sens
Sens.J                              =   []; % J - Jacobian
Sens.norm                           =   'r'; % either [], R, C or RC for non, row, column, row and column normalisation
Sens.C                              =   []; % C - column normalisation matrix (default Ikxk)
Sens.R                              =   []; % R; row normalisation matrix (default Imxm)

%==========================================================================%
% DEFINE INVERSE PROBLEM PARAMETERS
global Inv Sol
Inv.fwd_method                      =   'bicgstab'; % set forward solver method within inverse framework (either 'pcg', 'bicgstab', 'minres' or 'gmres' for pcg, stab BiCG, MinRes, GMRes)
Inv.fwd_tol                         =   1e-12;%1e-12; % forward solution within inverse framework relative residual tolerance (default 1e-12)
Inv.fwd_maxit                       =   10000;%100; % forward solution within inverse framework maximum number of iterations (default 100)
Inv.fwd_restarts                    =   20; % define restatrs for GMRes solver only (default 10)
Inv.fwd_precond                     =   'ilu'; % set preconditioner method (either 'chol', 'amg', 'ilu', 'camg' for incomple Cholesky, AMG, incomplete LU, complex AMG)
Inv.fwd_drop_tol                    =   1e-5; % drop tolerance fot incomplete preconditioners (default 1e-3)
Inv.linear                          =   true; % linear or non-linear inverse framework
Inv.method                          =   'tsvd'; % inverse method (either 'tsvd', 'tik', (if  Inv_lin is true) otherwise 'dgn', 'kndgn','lm', 'knlm', 'unlm', 'vmm', 'nlcg', 'cw' for TSVD, Tikhonov, Damped Guass-Newton, Krylov-Newton Damped Guass-Newton, Levenberg-Marqurdt, underdetermined Levenberg-Marqurdt, Krylov-Newton Levenberg-Marqurdt, Variable Metric, non-linear Conjugate Gradients, Compartment-wise reconstruction
Inv.GMRes_tol                       =   1e-6; % GMRes relative residual error tolerance for Krylov-Newton method only
Inv.TSVD_noise                      =   [];% [1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 ]; % TSVD noise criteria for determining which singular values attenuation will determine the truncation
Inv.maxit                           =   12; % maximum number of iterations for non-linear reconstruction
Inv.tol                             =   1e-5; % halting criteria for non-linear reconstruction  (defined by the differntial error)
Inv.regularisation_fn               =   []; % regularisation matrix file name
Inv.hyper_param                     =   'lc'; % regularisation hyper parameter selection approach (either 'lc', 'gcv', or 'kl' for L-curve, Generlised Cross Validation and KL distance)
Inv_grid                            =   []; % inversion grid dimensions (if empty no base transforamtion is performed)
Inv.regpar                          =   false; % selection of the regularisation parameter by Juanito: does Tikhonov inversion, GCV parameter selection, and Row normalisation. If false TSVD with fix truncation value will be used only instead
Inv.white                           =   false; % whitening of the noise

%==========================================================================%
% DEFINE SOLUTION PARAMETERS
Sol.ref                             =   mat_ref;%repmat(0.3,size(Mesh.tri,1),1); % S/m
Sol.current                         =   []; % sol; % current solution
Sol.rmask                           =   []; % rmask;  % enable focused reconstruction {ones for elements which requires reconstruction, 0 for the rest}
Sol.t_fac                           =   []; % tfac; % spatial regularisation hyper parameter
Sol.lamda                           =   []; % lamda; % trust region size
Sol.iter                            =   0; % iter; % iteration number
Sol.dur_it                          =   0; % dur_it % iteration duration
Sol.no_negs                         =   true; % f_neg; % flag to avoid negative conductivity values
Sol.rmask                           =   ones(size(Sol.ref));
Sol.rmask                           =   ones(size(Sol.ref));
rmask_ind                           =   find(Sol.rmask);

%==========================================================================%
% DEFINE DATA PARAMETERS
global Data
Data.fn                             =   '' ; % data file name
Data.source                         =   'Mk25'; % either 'Sim', 'Kr', 'Mk1b', 'Mk2' for Simulated, Korean system, Mark 1b, Mark 2/2.5
Data.pre                            =   false; % preprocessed data
Data.real                           =   true; % false
Data.freqs                          =   []; % frequencies by which the data was measured
Data.freqs.chosen                   =   1; % which frequencies to consider this is now taken from the view_data later
Data.frm.chosen                     =   []; % which frames to consider
Data.realibility                    =   []; % std / covariance inverse of the measured channels (set infinity for eliminated eletrodes) % could be over frequency as well
Data.comp_fn                        =   []; % compliance file name allow deriving contact impedance estimate
Data.bnd_v                          =   []; % =bnd_v; % no. of channels over time or frequency or both
Data.bnd_v_it                       =   []; % =bnd_v_it; % no. of channels over time or frequency or both
Data.description.r                  =   []; % =radii; % could be elipsoids as well
Data.description.perb_conductivity  =   []; % =sim_cnd;
Data.description.preb_location      =   []; % =location;
Data.description.preb_background    =   []; % =mat_ref;
Data.pca                            =   true; % PCA preprocessing by Juanito: it decomposes the data on the ppal components (PC) and projects the data of the first PCs, ie. removing noise.

%==========================================================================%
% DEFINE RASTERISATION PARAMETERS
global Image
Image.pixel_size                    =   3e-3;
Image.software                      =   'nim'; % 'graph', 'nim', 'maya' for rasterising for GraphEIT, Nim and MayaVI

%==========================================================================%
% DEFINE USER INTERACTION SETTINGS
global Progress
Progress.plots                      =   false; % print various plots during procedure
Progress.remarks                    =   false; % consule output regaring stages in which the code is
Progress.bars                       =   false;

% PROGRESS %
toc
display('SuperSolver starts...')
tic
%==========================================================================%
% MESH OPTIMISATION AND PROCESSING
[Mesh.vtx,Mesh.tri]                 =   optimize6_c(Mesh.vtx,Mesh.tri,Data.real,0);
[Mesh.srf]                          =   dubs3_2(Mesh.tri); % gets the elements on the outer surface of the mesh
[Mesh.sgrad, Mesh.volumes]          =   shape_functions_constructor(Mesh.vtx, Mesh.tri); % find the volumes of the elements

%==========================================================================%
% ASSIGN SURFACES TO EACH ELECTRODE %
[Fem.elec, Mesh.cnts, Fem.srf_area, Fem.center_of_mass]     =   set_electrodes_var_3(Mesh.vtx, Mesh.srf, Fem.pos, Fem.elec_diam, Fem.elec_method);

% print center of mass and area of electrode 30
fprintf('Electrode 30\nArea: %6.10f \nCenter of mass: %6.10f , %6.10f , %6.10f \n', Fem.srf_area(30), Fem.center_of_mass(30,1), Fem.center_of_mass(30,2), Fem.center_of_mass(30,3))


%==========================================================================%
% CHOOSE THE GROUND INDEX %
[Fem.gnd_ind]                       =   set_ground_force(Mesh.vtx, Mesh.srf, Fem.gnd_pos, Mesh.cnts,1);

%==========================================================================%
% SET CURRENTS %
[Fem.I, Fem.df]                     =   set_currents_2(Fem.prt, Fem.pos, Mesh.vtx, Fem.current);

%==========================================================================%
% SET CONTACT IMPEDANCE %
if length(Fem.zc) < 2
    [Fem.zc]                        =   contact_impedance(Data.comp_fn);
end

%==========================================================================%
% PROGRESS %
disp('Setup stage is completed')
disp('Forward modelling starts')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TIME IN ABOVE SECTION 16secs%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================%
% BOUNDARY VOLTAGE GENERATION %
fwd_parm_validator(); % validates that the parameters entered are valid

if SIM % CALCULATE FORWARD FOR PERTURBED CASE %

    Data.bnd_v_c0                   =   zeros(size(Fem.prt,1),BATCH_NO);
    Data.bnd_v_c                    =   zeros(size(Fem.prt,1),BATCH_NO);
    
    for batch_loop = 1:BATCH_NO
        tic
        Sol.ref_c                       =   1.1*Sol.ref;
        
        [Fem.E]                         =   fem_master_full_4(Sol.ref_c); % builds up system matrix based on complete electrode model
        Fwd.current_field               =   zeros(size(Fem.I)); % Initialize forward solution to zeros vector
        [Fwd]                           =   forward_solver_9(Fem.E,Fem.I); % Solve forward solution for the current field
        [Data.bnd_v_c(:,batch_loop)]    =   get_boundary_meas_2 (Fwd.current_field);
        Data.bnd_v_c0(:,batch_loop)     =   Data.bnd_v_c(:,batch_loop);
        toc
    end
    
    if ADD_NOISE % adds noise to the perturbation case
        Data.bnd_v_c  =   Data.bnd_v_c+Noise_Level*1e-6*randn(size(Data.bnd_v_c));
    end
    
end

toc


% set the perturbation

% strokelocation = [0.1262282,0.11792685,0.0958625];
% strokeradius = 10/1000.;
% Sol.ref_c = Sol.ref;
% count = 0;
% for i = 1:length(Sol.ref_c)
%     element_center = mean(Mesh.vtx(Mesh.tri(i,:),:));
%     if norm((element_center-strokelocation),2) < strokeradius
%         Sol.ref_c(i) = 0.9;%Sol.ref_c(i)*0.2;
%         count = count + 1;
%         %figure(99); hold on; plot3(element_center(1),element_center(2),element_center(3),'ro');
%     end
% end



% REPEAT FOR BASELINE/NON-PERTURBED CASE %
tic
[Fem.E]                             =   fem_master_full_4(Sol.ref);
%[Fem.E]                             =   fem_master_full_4(Sol.ref_c);
Fwd.current_field                   =   zeros(size(Fem.I));
toc
tic
[Fwd]                               =   forward_solver_9(Fem.E,Fem.I);
[Data.bnd_v]                        =   get_boundary_meas_2 (Fwd.current_field);
toc
tic
if SAVE
    save([Output_Data_Path tag], '-struct', 'Data', 'bnd*', '-v7.3')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TIME FOR EACH SIM 13 mins %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================%
% PROGRESS %
disp('finished forward')

if SOLVE % RUN THE REST OF THE FORWARD SOLVER AND INVERSION %
    
    %======================================================================%
    % SOLVE THE FORWARD SOLUTION FOR THE MEASUREMENT FIELD %
    Fwd.measurement_field           =   zeros(size(Fem.I,1),length(Fem.prt));
    [Fwd]                           =   adjoint_fields_9(Fem.E);
    toc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %TIME FOR FWD IS 52 mins %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %======================================================================%
    % CONSTRUCT THE JACOBIAN, THE FOCUSING IS DONE WHEN INVERTING STAGE %
    tic
    [Sens.J]                        =   jacobian_3d_7(Fwd.current_field, Fwd.measurement_field);%3.9mins
    toc
    qqqqq
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %TIME FOR SENS IS 8 mins %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %======================================================================%
    % PROGRESS %
    disp('finished J')
    
    %======================================================================%
    % INVERSION %
    for sv = 1:length(Inv.TSVD_noise)+1
        
        if sv>length(Inv.TSVD_noise) % fix no_svd to svd_list
            no_svd                  =   svd_list;
            Inv.TSVD_noise(sv)      =   0;
        else % variable no_svd (threshold method)
            disp(num2str([sv length(Inv.TSVD_noise)]))
            no_svd                  =   [];
        end
        
        tolerance                   =   0;
        normalise                   =   'y';%'y' for yes and '' for no
        f_support                   =   false;%true;
        
        if SIM
            diff_data0              =   100*( Data.bnd_v_c0 - Data.bnd_v)./Data.bnd_v;
            diff_data               =   100*( Data.bnd_v_c - Data.bnd_v)./Data.bnd_v;
        end
        
        if EXP_DATA
            diff_data_Exp           =   100*(Data.bnd_v_c_Exp - Data.bnd_v_Exp*ones(1,size(Data.bnd_v_c_Exp,2)))./(Data.bnd_v_Exp*ones(1,size(Data.bnd_v_c_Exp,2)));
        end
        
        for sv_loop=1:length(no_svd) % this not going to work anymore for rlative threshold -- just fixed
            
            tic
            [Jinv] =	TSVD_inverter_2 (Sens.J, Mesh.vtx, Mesh.tri, no_svd(sv_loop), tolerance, normalise, Sol.ref, Sol.rmask, Inv.TSVD_noise(sv), f_support, 0);
            time_Jinv=toc;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %TIME FOR JINV IS 51 mins%
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if SIM
                
                sol0                =   -Jinv*diff_data0; %noise free %5 mins
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %TIME FOR JINV IS 4 mins %
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if ADD_NOISE
                    
                    sol             =   -Jinv*diff_data;% with noise
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    %TIME FOR JINV IS 8 mins %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                end
            end
            
            if EXP_DATA
                sol_Exp             =   -Jinv*diff_data_Exp;%(:,1:DiluteFactor:end);% experimental data
            end
            
            if SAVE
                
                if exist('sol','var')
                    save([Output_Sol_Path '_' num2str(no_svd(sv_loop)) 'SVs'], 'sol*','Jinv','Data','-v7.3');
                else
                    save([Output_Sol_Path '_' num2str(no_svd(sv_loop)) 'SVs'], 'Jinv','Data','-v7.3');
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %TIME FOR JINV IS 7 mins %
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
            
            clear sol* Jinv
            
        end % END OF 'SV_LOOP' LOOP %
        
    end % END OF 'SV' LOOP %
    
end % END OF IF SOLVE %

%clear all;close all;clc