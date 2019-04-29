function fwd_parm_validator(Fwd,Fem)
% Usage: fwd_parm_validator();
%
% General: validates the forward parameters
%
% Input:
%
% Output:
%
% LH 15/12/05
%==============================================================================%
% global Progress Fwd Fem

switch lower(Fwd.method)
    case {'stabbicg','minres','gmres'}
        if ~strcmp(lower(Fwd.precond),'ilu')
            disp_progress('Preconditioning method is not available yet, setting preconditioner to ILU')
            Fwd.precond = 'ilu';
        end
end
% set minimum residual error tolerance to 1e-16
if Fwd.tol < 1e-16
    disp_progress ('Forward model tolerance is too low setting to 1e-16')
    Fwd.tol = 1e-16;
end
