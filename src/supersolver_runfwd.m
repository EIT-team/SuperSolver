function [Mesh,Fem,Fwd,Inv,Sol,Data] = supersolver_runfwd(Mesh,Fem,Fwd,Inv,Sol)
%supersolver_runfwd Summary of this function goes here
%   Detailed explanation goes here


disp('Supersolver running...');

fwd_parm_validator(Fwd,Fem); % validates that the parameters entered 

[E,Fem,Mesh]                             =   fem_master_full_4(Sol.ref,Fem,Mesh);
Fem.E=E;
Fwd.current_field                   =   zeros(size(Fem.I));
[Fwd]                               =   forward_solver_9(Fwd,Fem);
[Data.bnd_v]                        =   get_boundary_meas_2 (Fwd,Fem,Mesh);

end

