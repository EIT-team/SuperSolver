%% Initialise
disp('Setup started')
tic
% load existing mesh
M=load('..\resources\SA060.mat'); %cyl tank 32 chn from zz
%create conductivity vector
mat_ref=M.mat_ref;
mat_ref(mat_ref==1) = 0.35;
% NEED TO VALIDATE MAT_REF = N x1

% load protocol
P=load('..\resources\SA060_prt.mat');

[Mesh,Fem,Fwd,Inv,Sol] = supersolver_init(M.vtx,M.tri,mat_ref,M.pos,M.gnd_pos,P.prt_full);

toc
disp('SuperSolver starts...')

%% Setup

[Mesh,Fem,Fwd,Inv,Sol] = supersolver_setup(Mesh,Fem,Fwd,Inv,Sol);

toc
disp('Setup stage is completed')
disp('Forward modelling starts')

%% run fwd


[Mesh,Fem,Fwd,Inv,Sol,Data] = supersolver_runfwd(Mesh,Fem,Fwd,Inv,Sol);
toc
disp('finished forward')

%% get jac

[Mesh,Fem,Fwd,Inv,Sol,J] = supersolver_solve(Mesh,Fem,Fwd,Inv,Sol);

toc
disp('finished J')

