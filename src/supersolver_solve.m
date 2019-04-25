function [Mesh,Fem,Fwd,Inv,Sol,J] = supersolver_solve(Mesh,Fem,Fwd,Inv,Sol)
%supersolver_solve Creates the jacobian for given forward solution
%   Detailed explanation goes here

%======================================================================%
% SOLVE THE FORWARD SOLUTION FOR THE MEASUREMENT FIELD %
Fwd.measurement_field           =   zeros(size(Fem.I,1),length(Fem.prt));
[Fem,Fwd,Inv]                   =   adjoint_fields_9(Fem,Fwd,Inv,Fem.E);

[J]                        =   jacobian_3d_7(Fem,Fwd,Mesh,Fwd.current_field, Fwd.measurement_field);

end

