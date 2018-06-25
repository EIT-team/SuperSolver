%clear all;close all;%clc
clear Mesh;
disp('Setup started')
tic
%==========================================================================%
% SWITCH SETUP SECTION %
ADD_NOISE                   =       0;
Noise_Level_additive        =       2; % noise level in uV
Noise_Level_proportional    =    0.02; % proportional noise level in %
SIM                         =       0;  % 0 - no pert, 1 - with pert
FUCK_YOU                    =       0;

SOLVE       =       0;  % 0 - stops after boundary voltage generation, 1 - runs through to inversion
PLOT        =       0;

%==========================================================================%
% DEFINE MESH PARAMETERS %
global Mesh
Mesh.fn                             =   'Cylindrical Tank';
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


% load existing mesh
M=load('..\resources\SA0602.mat'); %cyl tank 32 chn from zz

% adjust geometry here


%set geometry
Mesh.vtx = M.vtx;
Mesh.tri = M.tri;

%set conductivities
mat_ref(M.mat_ref ==1) = 0.35;

%==========================================================================%
% DEFINE FEM PARAMETERS
global Fem
Fem.current                         =   300e-6; % 50e-6;% [50uA in rat and tank expt, 400uA is some simulations]
Fem.elec_diam                       =   11.5e-3;% electrode diameter in meters
Fem.materials_fn                    =   [];
Fem.elec_method                     =   's'; % electrodes assigning method 'n', 's', 'm', 'v' closest node, closest elemt centre, nodes within Elec_diam/2 from remote-most surface coordinate, nodes within Elec_diam/2
Fem.df                              =   []; % df;
Fem.elec                            =   []; % elec;
Fem.gnd                             =   []; % gnd_ind;
Fem.I                               =   []; % I - sources matrix
Fem.zc                              =   200 ; % contact impedance (could be a vector at the length of the number of electrodes);
Fem.E                               =   []; % E; FEM system matrix
Fem.isotropy                        =   true; % define whether a model is Isotropic or anisotropic

% ==========================================================================%
% DEFINE THE PROTOCOL %
P=load('..\resources\SA060_prt.mat');

%adjust positions here

Fem.prt = P.prt_full;

Fem.pos                             =   M.pos;
Mesh.pos                            =   M.pos;
Fem.gnd_pos                         =   M.gnd_pos; %SA060 GND


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


toc
disp('SuperSolver starts...')
tic
%==========================================================================%
% MESH OPTIMISATION AND PROCESSING
[Mesh.vtx,Mesh.tri]                 =   optimize6_c(Mesh.vtx,Mesh.tri,true,0);
[Mesh.srf]                          =   dubs3_2(Mesh.tri); % gets the elements on the outer surface of the mesh
[Mesh.sgrad, Mesh.volumes]          =   shape_functions_constructor(Mesh.vtx, Mesh.tri); % find the volumes of the elements

%==========================================================================%
% ASSIGN SURFACES TO EACH ELECTRODE %
[Fem.elec, Mesh.cnts, Fem.srf_area, Fem.center_of_mass]     =   set_electrodes_var_3(Mesh.vtx, Mesh.srf, Fem.pos, Fem.elec_diam, Fem.elec_method);

%==========================================================================%
% CHOOSE THE GROUND INDEX %
[Fem.gnd_ind]                       =   set_ground_force(Mesh.vtx, Mesh.srf, Fem.gnd_pos, Mesh.cnts,1);

%==========================================================================%
% SET CURRENTS %
[Fem.I, Fem.df]                     =   set_currents_2(Fem.prt, Fem.pos, Mesh.vtx, Fem.current);

%==========================================================================%
% SET CONTACT IMPEDANCE %
if length(Fem.zc) < 2
    [Fem.zc]                        =   contact_impedance([]);
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

if SIM % CALCULATE FORWARD FOR PERTURBED CASE
    %% set the perturbation
    
    strokelocation = [0.13, 0.16, 0.11];
    %strokelocation = [0.1262282,0.11792685,0.0958625];
    %strokelocation = [0,0,0];
    strokeradius = 20/1000.;
    Sol.ref_c = Sol.ref;
    count = 0;
    for i = 1:length(Sol.ref_c)
        element_center = mean(Mesh.vtx(Mesh.tri(i,:),:));
        if norm((element_center-strokelocation),2) < strokeradius
            Sol.ref_c(i) = 0.7;%Sol.ref_c(i)*0.2;
            count = count + 1;
            figure(1); hold on; plot3(element_center(1),element_center(2),element_center(3),'ro');
        end
    end
    
    tic
    [Fem.E]                         =   fem_master_full_4(Sol.ref_c); % builds up system matrix based on complete electrode model
    Fwd.current_field               =   zeros(size(Fem.I)); % Initialize forward solution to zeros vector
    [Fwd]                           =   forward_solver_9(Fem.E,Fem.I); % Solve forward solution for the current field
    Data.bnd_v_c    =   get_boundary_meas_2 (Fwd.current_field);
    Data.bnd_v_c0     =   Data.bnd_v_c;
    toc
    
    if ADD_NOISE % adds noise to the perturbation case
        Data.bnd_v_c  =   (1+randn(size(Data.bnd_v_c))*Noise_Level_proportional/100).*Data.bnd_v_c+Noise_Level_additive*1e-6*randn(size(Data.bnd_v_c));
    end
    
end

toc


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

disp('finished forward')

%% START RECONSTRUCTION

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %TIME FOR SENS IS 8 mins %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %======================================================================%
    % PROGRESS %
    disp('finished J')
    
    %====== interrupt here if no reconstruction should be made ======%
    %stop
    %================================================================%
    
    diff_data0 = Data.bnd_v_c0 - Data.bnd_v;
    diff_data = Data.bnd_v_c - Data.bnd_v;
    
    %% ZERO ORDER TIKHONOV
    
    load('.\mesh_prot_elec\TA052_hex.mat');
    J_hex = jacobian_hexagoniser(Sens.J,TA052_hex);
    
    lambda = logspace(-12,0,3000);
    
    [U,sm,X,V,W] = cgsvd(J_hex,eye(size(J_hex,2),size(J_hex,2)));
    [x_lambda_0,rho_0,eta_0] = tikhonov(U,sm,X,diff_data,lambda);
    
    % find L-curve corner and thus the optimal lambda
    [reg_c,rho_c,eta_c] = l_corner(rho_0,eta_0,lambda);
    [min_lambda0,min_index0] = min(abs(lambda-ones(size(lambda))*reg_c));
    
    disp('done with zero order tikhonov');
    
    %% FIRST ORDER TIKHONOV
    disp('starting first order tikhonov');
    
    xyz=(TA052_hex.Nodes(TA052_hex.Hex(:,1),:)+TA052_hex.Nodes(TA052_hex.Hex(:,2),:)+TA052_hex.Nodes(TA052_hex.Hex(:,3),:)+...
        TA052_hex.Nodes(TA052_hex.Hex(:,4),:)+TA052_hex.Nodes(TA052_hex.Hex(:,5),:)+TA052_hex.Nodes(TA052_hex.Hex(:,6),:))./6;
    L = gen_Laplacian_matrix(xyz,1e-2);
    [U,sm,X,V,W] = cgsvd(J_hex,L);
    [x_lambda_1,rho_1,eta_1] = tikhonov(U,sm,X,diff_data,lambda);
    
    % find L-curve corner and thus the optimal lambda
    [reg_c,rho_c,eta_c] = l_corner(rho_1,eta_1,lambda);
    [min_lambda1,min_index1] = min(abs(lambda-ones(size(lambda))*reg_c));
    
    disp('done with first order tikhonov');
    
    
end