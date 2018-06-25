function [Fwd, flag, relres, iter] = adjoint_fields_9(E, current, tol, method);
%Usage: [Fwd, flag, relres, iter] = adjoint_fields_9(E, current,tol, method);
%
% General:
% real conductivities - AMG preconditioner (based on Powell CE (FEMLAB) code) and PCG solver  or Cholesky (peconditioner) and PCG
% complex admittivities - LU (preconditioner) stabilized Bi-CG or GMRES
% ver 6. allow setting a specific measurement "current" and choosing the desired solver
% ver. 8 % ver. 8 allow globals to define preconditioning and solution method
% Important Note: this function assumes that v_f and dur_it have already externally defined as global of the following dimentions
% {n + no. of electrodes x no. of measurements} and a {scalar} respectivly
% ver. 9 assumes globals
%
% Note:
% Matlab R14 7.0.1, has a bug in gmres.m function. in case of an error
% message, please follow the instructions at:
% http://www.mathworks.com/support/solutions/data/1-Z81WY.html?solution=1-Z81WY
%
%
% Input:
% E - full rank system matrix based on the 3D complete electrode model {n + no. of electrodes x n + no. of electrodes}
% current - injected "measurement" current (in Amps)
% m_ind (global) - measurements matrix (indices of electrode pairs)
%
% v_f (global) -  approximated nodal potential distribution of the adjoint field  {n + no. of electrodes x no. of measurements}
% Inv.fwd_method (global) - set forward solver method (either 'pcg','bicgstab','minres' or 'gmres' for PCG, stab BiCG, MinRes, GMRes)
% Inv.fwd_tol (global) - fwd solution relative residual tolerance (default 1e-12) {scalar}
% Inv.fwd_drop_tol (global) - for incomplete factorisation {scalar}
% Inv.precond (global) - set preconditioner method (either 'chol', 'amg', 'ilu', 'camg' for incomple Cholesky, AMG, incomplete LU, complex AMG)
% Inv.fwd_maxit (global) - maximum iterations for the forward solver {scalar}
% Inv.restarts (global)- defines after how many vectors GMRes retains before it restarts itself
%
% tol (for backwards support) - forward solution tolerance,  or alternatively includes preconditioner drop tolerance as a second entry {scalar} / {2 x 1}
% fmethod (for backwards support) - allow the following
%          'a' - AMG preconditioner and PCG
%          'c' - Cholesky and PCG
%          'b' - LU and stabilized Bi-CG
%          'g' - LU and GMRES
%
% Output:
% v_f - measurement field {n + no. of electrodes x no. of measurements}
% flag - convergence flags vector {number of independent current patterns x 1}
%    0 converged to the desired tolerance tol within maxit iterations
%    1 iterated maxit times but did not converge
%    2 preconditioner was ill-conditioned
%    3 stagnated solution (two consecutive iterates were the same)
% relres - relative residual error (i.e. the norm ||Ax-b||/||b||) {number of independent current patterns x 1}
% iter - outer and inner (for gmres) iterations {number of independent current patterns x 1 (or 2 for gmres)}
%
% LH 15/12/05
%------------------------------------------------------------------------------------------------------------------------
% load the previous solution
global Fem Fwd Inv

% in case v_f was not defined externally
if isempty(Fwd.measurement_field)
    Fwd.measurement_field  = zeros(length(Fem.E),length(Fem.prt));
    disp_progress('Fwd.measurement_field was not defined as global before setting it to zero')
end

% to be consistent with older versions
if nargin < 2
    current = 1; % 1 Amp   
end

% in case the user would like to override the chosen setting by the globals
if nargin > 3
    if length(tol) < 2
        Inv.fwd_drop_tol = 1e-5;% set preconditioner default tolerance in case it was not indicated
    end
    Inv.fwd_tol = tol(1);
    Inv.fwd.restarts = 20;
    Inv.fwd_maxit = 100;
    switch lower(method)
        case {'amg','a'}
            Inv.fwd_precond = 'amg';
            Inv.fwd_method = 'pcg';
        case {'pcg','c'}
            Inv.fwd_precond = 'chol';
            Inv.fwd_method = 'pcg';
        case {'bicg','b'}
            Inv.fwd_precond = 'ilu';
            Inv.fwd_method = 'bicgstab';
        otherwise
            Inv.fwd_precond = 'ilu';
            Inv.fwd_method = 'gmres'
    end
else % accept the globals settings
     E = Fem.E;
end

% this bit can be calculated only once externally
[unique_m_ind, un_i, un_j] = unique(Fem.prt(:,3:4), 'rows'); % pick only the unique measurement patterns
Is_supl = zeros(length(E) - size(Fem.elec,1), size(unique_m_ind,1));
%no of electrodes x no of measurements (now currents)!

MC = [];
for ii=1:size(unique_m_ind,1)
    m_n = zeros(size(Fem.elec,1),1);
    m_n(unique_m_ind(ii,1)) = current;
    m_n(unique_m_ind(ii,2)) = -current;
    MC = [MC, m_n];
end

I = [Is_supl; MC];
I(Fem.gnd_ind,:) = 0;

unique_v_f = zeros(size(E,1),size(MC,2));

% tic
switch lower(Inv.fwd_method)
    case {'pcg','p'}
        if strcmp(lower(Inv.fwd_precond),'amg')
            options.mi=2;
            options.spc=1;
            [m1,m2] = amgset(E, options); % set the grids and preconditioning
            %         dur_precond = toc;
            for ii=1:size(MC, 2)
                [unique_v_f(:,ii), flagp(ii), relres(ii), iter(ii), resvec] = pcg(E,I(:,ii),Inv.fwd_tol,Inv.fwd_maxit,@amgsol,[],Fwd.measurement_field  (:,un_i(ii)),m2);
            end
            %         dur_po = toc;
            Fwd.measurement_field   = unique_v_f(:,un_j);
            %         dur_it = dur_it + dur_po;
            % un-comment this part in case you would like to track the accuracy and duration of the forward solution through the iterations
            % res_po = max(relres(ii));
            % disp(['m3d AMG - relres ', num2str(res_po),'   duration ',num2str(dur_po), ' + ',num2str(dur_precond) ,'   Inv.fwd_tol ', num2str(Inv.fwd_tol),'   iter ', num2str(mean(iter))])
        elseif strcmp(lower(Inv.fwd_precond),'chol')
            %disp('Method is Cholesky and PCG')
           K = cholinc(E, Inv.fwd_drop_tol); % sparse Cholesky factorization
            %         dur_precond = toc;
            for ii=1:size(MC, 2)
                [unique_v_f(:,ii), flagp(ii), relres(ii), iter(ii), resvec] = pcg(E,I(:,ii),Inv.fwd_tol,Inv.fwd_maxit,K',K,Fwd.measurement_field  (:,un_i(ii)));
            end
            %         dur_po = toc;
            Fwd.measurement_field   = unique_v_f(:,un_j);
            %         dur_it = dur_it + dur_po;
            % un-comment this part in case you would like to track the accuracy and duration of the forward solution through the iterations
            % res_po = max(relres(ii));
            % disp(['m3d Cholesky - relres ', num2str(res_po),'   duration ',num2str(dur_po), ' + ',num2str(dur_precond) ,'   Inv.fwd_tol ', num2str(Inv.fwd_tol),'   iter ', num2str(mean(iter))])
        end
        return
    case {'bicgstab','b'}
        %disp('Method is LU and Stabilized -BiCG')
        %[L,U] = luinc(E, Inv.fwd_drop_tol); % sparse LU decompositon
        %         dur_precond = toc;
 %       for ii=1:size(MC, 2)
            %[unique_v_f(:,ii),flag(ii),relres,iter,resvec] = bicgstab(E,I(:,ii),Inv.fwd_tol,Inv.fwd_maxit,L,U,Fwd.measurement_field  (:,un_i(ii)));
            %[unique_v_f(:,ii),flag(ii),relres,iter,resvec] = bicgstab(E,I(:,ii),Inv.fwd_tol,Inv.fwd_maxit,diag(diag(E)));
 %           [Fwd.current_field(:,ii),flag(ii),relres(ii),iter(ii),resvec] = pcg(E,I(:,ii), Inv.fwd_tol,Inv.fwd_maxit,diag(diag(E)));
 %       end
        %         dur_po = toc;
 %       Fwd.measurement_field   = unique_v_f(:,un_j);
        %         dur_it = dur_it + dur_po;
        % un-comment this part in case you would like to track the accuracy and duration of the forward solution through the iterations
 %       res_po = max(relres(ii));
        
        %markus
        setup.droptol = Inv.fwd_drop_tol;
        [L,U] = ilu(E,setup);
        disp('start bicgstab');
        for ii=1:size(MC,2)
            [unique_v_f(:,ii),flag(ii),relres(ii),iter(ii),resvec]=pcg(E,I(:,ii),Inv.fwd_tol,Inv.fwd_maxit,L,U,Fwd.measurement_field(:,un_i(ii)));
            if flag(ii)
                disp('pcg did not converge!!');
            end
        end
        
        Fwd.measurement_field = unique_v_f(:,un_j);

        return
    case {'minres','m'}
        %disp('Method is LU and Stabilized -BiCG')
        [L,U] = luinc(E, Inv.fwd_drop_tol); % sparse LU decompositon
        %         dur_precond = toc;
        for ii=1:size(MC, 2)
            [unique_v_f(:,ii),flag(ii),relres,iter,resvec] = minres(E,I(:,ii),Inv.fwd_tol,Inv.fwd_maxit,L,U,Fwd.measurement_field  (:,un_i(ii)));
        end
        %         dur_po = toc;
        Fwd.measurement_field   = unique_v_f(:,un_j);
        %         dur_it = dur_it + dur_po;
        % un-comment this part in case you would like to track the accuracy and duration of the forward solution through the iterations
        % res_po = max(relres(ii));
        % disp(['m3d LU+BiCG - relres ', num2str(res_po),'   duration ',num2str(dur_po), ' + ',num2str(dur_precond) ,'   Inv.fwd_tol ', num2str(Inv.fwd_tol),'   iter ', num2str(mean(iter))])
        return
    otherwise
        % disp('Method is GMRES')
        [L,U] = luinc(E, Inv.fwd_drop_tol);% sparse LU decompositon
        %         dur_precond = toc;
        for ii=1:size(MC, 2)
            [unique_v_f(:,ii),flag(ii),relres,iter,resvec] = gmres(E,I(:,ii), Inv.fwd_restarts,Inv.fwd_tol,Inv.fwd_maxit,L,U,Fwd.measurement_field  (:,un_i(ii)));
        end
        %         dur_po = toc;
        Fwd.measurement_field   = unique_v_f(:,un_j);
        %         dur_it = dur_it + dur_po;
        % un-comment this part in case you would like to track the accuracy and duration of the forward solution through the iterations
        % res_po = max(relres(ii));
        % disp(['m3d LU+GMRes - relres ', num2str(res_po),'   duration ',num2str(dur_po), ' + ',num2str(dur_precond) ,'   Inv.fwd_tol ', num2str(Inv.fwd_tol),'   iter ', num2str(mean(iter))])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This might become a part of the EIDORS suite
% Copyright (c) Lior Horesh 2005, EIT group, Medical Physics and Bioengineering, UCL, UK
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details
% EIDORS 3D version 2
% MATLAB Version 6.1.0.450 (R12.1) on PCWIN
% MATLAB License Number: 111672
% Operating System: Microsoft Windows Server 2003 Standard Edition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%