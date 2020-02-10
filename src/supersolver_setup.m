function [Mesh,Fem,Fwd,Inv,Sol] = supersolver_setup(Mesh,Fem,Fwd,Inv,Sol,PlotFlag)
%supersolver_setup Setup boundary conditions system matrix and sets up supersolver
%   Detailed explanation goes here

if nargin <6
    PlotFlag=0;
end


disp('Supersolver setup...'); 
%==========================================================================%
% MESH OPTIMISATION AND PROCESSING
[Mesh.vtx,Mesh.tri]                 =   optimize6_c(Mesh.vtx,Mesh.tri,true,PlotFlag);
[Mesh.srf]                          =   dubs3_2(Mesh.tri); % gets the elements on the outer surface of the mesh
[Mesh.sgrad, Mesh.volumes]          =   shape_functions_constructor(Mesh.vtx, Mesh.tri); % find the volumes of the elements

%==========================================================================%
% ASSIGN SURFACES TO EACH ELECTRODE %
[Fem.elec, Mesh.cnts, Fem.srf_area, Fem.center_of_mass]     =   set_electrodes_var_3(Mesh.vtx, Mesh.srf, Fem.pos, Fem.elec_diam, Fem.elec_method,PlotFlag);

%==========================================================================%
% CHOOSE THE GROUND INDEX %
[Fem.gnd_ind]                       =   set_ground_force(Mesh.vtx, Mesh.srf, Fem.gnd_pos, Mesh.cnts,PlotFlag);

%==========================================================================%
% SET CURRENTS %
[Fem.I, Fem.df]                     =   set_currents_2(Fem.prt, Fem.pos, Mesh.vtx, Fem.current);

%==========================================================================%
% SET CONTACT IMPEDANCE %
if length(Fem.zc) < 2
%     [Fem.zc]                        =   contact_impedance(Fem);
    Fem.zc = Fem.zc * ones(size(Fem.pos,1),1);
end


end

