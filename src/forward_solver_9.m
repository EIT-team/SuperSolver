function [Fwd, flag, relres, iter] = forward_solver_9(Fwd,Fem, tol, fmethod)
% Usage: [Fwd, flag, relres, iter] = forward_solver_9(E, I, tol, fmethod);
%
% General:
% This function solves the forward problem using the following methods:
% real conductivities - AMG preconditioner (based on Powell CE (FEMLAB) code) and PCG solver or Cholesky (peconditioner) and PCG
% complex admittivities - LU (preconditioner) stabilized Bi-CG or GMRES
%
% ver. 7c output conversgence parameters
% ver. 8 allow globals to define preconditioning and solution method
% ver. 9 assumes globals
%
% Note:
% Matlab R14 7.0.1, has a bug in gmres.m function. in case of an error
% message, please follow the instructions at:
% http://www.mathworks.com/support/solutions/data/1-Z81WY.html?solution=1-Z81WY
%
% Input:
% E - full rank system matrix based on the 3D complete electrode model {n + no. of electrodes x n + no. of electrodes}
% I -  currents matrix (RHS)  {number of independent current patterns x n + no. of electrodes}
%
% V (global) -  approximated nodal potential distribution {number of independent current patterns x n + no. of electrodes}
% Fwd.method (global) - set forward solver method (either 'pcg','bicgstab','minres' or 'gmres' for PCG, stab BiCG, MinRes, GMRes)
% Fwd.tol (global) - fwd solution relative residual tolerance (default 1e-12) {scalar}
% Fwd.drop_tol (global) - for incomplete factorisation {scalar}
% Fwd.precond (global) - set preconditioner method (either 'chol', 'amg', 'ilu', 'camg' for incomple Cholesky, AMG, incomplete LU, complex AMG)
% Fwd.maxit (global) - maximum iterations for the forward solver {scalar}
% Fwd.restarts (global)- defines after how many vectors GMRes retains before it restarts itself
%
% tol (for backwards support) - forward solution tolerance,  or alternatively includes preconditioner drop tolerance as a second entry {scalar} / {2 x 1}
% fmethod (for backwards support) - allow the following
%          'a' - AMG preconditioner and PCG
%          'c' - Cholesky and PCG
%          'b' - LU and stabilized Bi-CG
%          'g' - LU and GMRES
%
% Output:
% V (as global) -  approximated nodal potential distribution
% flag - convergence flags vector {number of independent current patterns x 1}
%    0 converged to the desired tolerance tol within maxit iterations
%    1 iterated maxit times but did not converge
%    2 preconditioner was ill-conditioned
%    3 stagnated solution (two consecutive iterates were the same)
% relres - relative residual error (i.e. the norm ||Ax-b||/||b||) {number of independent current patterns x 1}
% iter - outer and inner (for gmres) iterations {number of independent current patterns x 1 (or 2 for gmres)}
%
% LH 10/2/04, last updated 15/12/05
%------------------------------------------------------------------------------------------------------------------------
% load previous solution

% in case the user would like to override the chosen setting by the globals
if nargin > 2
    if length(tol)<2
        Fwd.drop_tol = 1e-5;% set preconditioner default tolerance in case it was not indicated
    end
    Fwd.tol = tol(1);
    Fwd.restarts = 20;
    Fwd.maxit = 100;
    switch lower(method)
        case {'amg','a'}
            Fwd.recond = 'amg';
            Fwd.method = 'pcg';
        case {'pcg','c'}
            Fwd.precond = 'chol';
            Fwd.method = 'pcg';
        case {'bicg','b'}
            Fwd.precond = 'ilu';
            Fwd.method = 'bicgstab';
        otherwise
            Fwd.precond = 'ilu';
            Fwd.method = 'gmres';
    end
else % accept the globals settings
    E = Fem.E;
    I = Fem.I;
end

d = size(I,2);% number of unique current patterns

% in case V was not defined externally
if isempty(Fwd.current_field)
    Fwd.current_field = zeros(size(E,1),d);
    disp_progress('Fwd.current_field was not defined as global before setting it to zero')
end

% tic
switch lower(Fwd.method)
    case {'pcg','p'}
        if strcmp(lower(Fwd.precond),'amg')
            %disp('Method is AMG and PCG')
            options.mi=2;
            options.spc=1;

            [m1,m2] = amgset(E, options); % set the grids and preconditioning
            for ii=1:d
                [Fwd.current_field(:,ii),flag(ii),relres(ii),iter(ii),resvec] = pcg(E,I(:,ii), Fwd.tol,Fwd.maxit,@amgsol,[],Fwd.current_field(:,ii),m2);
            end
            %         dur_po = toc;
            %         dur_it = dur_it + dur_po;
        elseif strcmp(lower(Fwd.precond),'chol')
            %disp('Method is Cholesky and PCG')
            K = cholinc(E, Fwd.drop_tol); % sparse Cholesky factorization
            for ii=1:d
                [Fwd.current_field(:,ii),flag(ii),relres(ii),iter(ii),resvec] = pcg(E,I(:,ii), Fwd.tol,Fwd.maxit,K',K,Fwd.current_field(:,ii));
            end
            %         dur_po = toc;
            %         dur_it = dur_it + dur_po;
        end
        return
    case {'bicgstab','b'}
        %disp('Method is LU and Stabilized -BiCG')
%         tic
%         [L,U] = luinc(E, Fwd.drop_tol); % sparse LU decompositon
%         toc
%         disp('start bicgstab');
%         for ii=1:d
%             tic
%             [Fwd.current_field(:,ii),flag(ii),relres(ii),iter(ii),resvec]=pcg(E,I(:,ii),Fwd.tol,Fwd.maxit,L,U,Fwd.current_field(:,ii)); %bicgstab
%             toc
%         end
        setup.droptol = Fwd.drop_tol;
        [L,U] = ilu(E,setup);
        disp('start pcg');
        for ii=1:d
            [Fwd.current_field(:,ii),flag(ii),relres(ii),iter(ii),resvec]=pcg(E,I(:,ii),Fwd.tol,Fwd.maxit,L,U,Fwd.current_field(:,ii));
            if flag(ii)
                disp('pcg did not converge!!');
            end
        end
        %         dur_po = toc;
        %         dur_it = dur_it + dur_po;
        return
    case {'minres','m'}
        %disp('Method is LU and MinRes')
        [L,U] = luinc(E, Fwd.drop_tol); % sparse LU decompositon
        for ii=1:d
            [Fwd.current_field(:,ii),flag(ii),relres(ii),iter(ii),resvec] = minres(E,I(:,ii), Fwd.tol,Fwd.maxit,L,U,Fwd.current_field(:,ii));
        end
        %         dur_po = toc;
        %         dur_it = dur_it + dur_po;
        return
    otherwise
        % disp('Method is GMRES')
        [L,U] = luinc(E, Fwd.drop_tol);% sparse LU decompositon
        for ii=1:d
            [Fwd.current_field(:,ii),flag(ii) ,relres(ii,:), iter(ii,:), resvec] = gmres(E,I(:,ii), Fwd.restarts, Fwd.tol, Fwd.maxit,L,U,Fwd.current_field(:,ii));
        end
        %         dur_po = toc;
        %         dur_it = dur_it + dur_po;
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