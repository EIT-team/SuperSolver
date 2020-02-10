function [Mesh,Fem,Fwd,Inv,Sol] = supersolver_init(vtx,tri,mat_ref,pos,gnd_pos,prt)
%supersolver_init Creates the structures needed for Supersolver 
%   Detailed explanation goes here

disp('Supersolver initialising...');

%% Create Structures
% some of these are overwritten but its easier having them all in one place

%==========================================================================%
% DEFINE MESH PARAMETERS %%
Mesh.fn                             =   'Generic Mesh';
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
Mesh.pos                            =   []; % electrode positions in same units as mesh

%==========================================================================%
% DEFINE FEM PARAMETERS
Fem.current                         =   300e-6; %injected current in uA
Fem.elec_diam                       =   10e-3;% electrode diameter in meters
Fem.materials_fn                    =   [];
Fem.elec_method                     =   's'; % electrodes assigning method 'n', 's', 'm', 'v' closest node, closest elemt centre, nodes within Elec_diam/2 from remote-most surface coordinate, nodes within Elec_diam/2
Fem.df                              =   []; % df;
Fem.elec                            =   []; % elec;
Fem.gnd                             =   []; % gnd_ind;
Fem.I                               =   []; % I - sources matrix
Fem.zc                              =   200 ; % contact impedance (could be a vector at the length of the number of electrodes);
Fem.E                               =   []; % E; FEM system matrix
Fem.isotropy                        =   true; % define whether a model is Isotropic or anisotropic
Fem.prt                             =   []; % measurement protocol I+ I- V+ V-
Fem.pos                             =   []; % electrode positions in same units as mesh
Fem.gnd_pos                         =   []; % postion of ground node NOT THE GROUND ELECTRODE

%==========================================================================%
% DEFINE FORWARD SOLVER PARAMETERS
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
% Sens.J                              =   []; % J - Jacobian
% Sens.norm                           =   'r'; % either [], R, C or RC for non, row, column, row and column normalisation
% Sens.C                              =   []; % C - column normalisation matrix (default Ikxk)
% Sens.R                              =   []; % R; row normalisation matrix (default Imxm)

%==========================================================================%
% DEFINE INVERSE PROBLEM PARAMETERS
% global Inv Sol
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
Inv.grid                            =   []; % inversion grid dimensions (if empty no base transforamtion is performed)
Inv.regpar                          =   false; % selection of the regularisation parameter by Juanito: does Tikhonov inversion, GCV parameter selection, and Row normalisation. If false TSVD with fix truncation value will be used only instead
Inv.white                           =   false; % whitening of the noise

%==========================================================================%
% DEFINE SOLUTION PARAMETERS
Sol.ref                             =   repmat(0.3,size(Mesh.tri,1),1); % S/m - generic conductivity to start with
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

%% Adjust Geometry
%set geometry
Mesh.vtx = vtx;
Mesh.tri = tri;

%set conductivities
% NEED TO VALIDATE MAT_REF = N x1
Sol.ref                             =   double(mat_ref); % S/m ensure double

%% Adjust Electrodes and Protocol

Fem.pos   =   pos;
Mesh.pos  =   pos;
Fem.gnd_pos =   gnd_pos;

Fem.prt = prt;

end

