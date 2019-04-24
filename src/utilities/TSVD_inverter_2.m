function [Jinv] = TSVD_inverter_2 (J, vtx, tri, no_svd, tolerance, normalise, mat_ref, rmask, noise, f_support, PLOT)
%function [Jinv,no_svd] = TSVD_inverter_2 (J, vtx, tri, no_svd, tolerance, normalise, mat_ref, rmask, noise, f_support)
% Usage: [Jinv] = TSVD_inverter_2 (J, vtx, tri, no_svd, tolerance, normalise, mat_ref, rmask, noise, f_support);
%
% General:
% performs TSVD pseudo-inversion
% applicable for both nodal and element basis
% ver. 2 allow discarding support from reconstruction
%
% Input:
% J - Jacobian {m x k} or {m x n}
% vtx - vertices {n x 3}
% tri - simplices {k x 4}
% no_svd - truncation point for SVD (if provided) {scalar}
% tolerance - truncation tolerance {scalar}
% normalise - flag for row and column normalisation {boolean}
% mat_ref - admittivity refernce {k x 1}
% rmask - ROI mask {k x 1} or {n x 1}
% noise - trunction ratio level {scalar}
% f_support - flag for applying support {boolean}
% PLOT - flag for plotting, 0 - is off, 1 - is on
% 
% Output:
% Jinv - sC(RAC)^ penrode moore pseudo inverted Jacobian
%%% compensate - ratio of SV to truncated SV {scalar}
% based on old Liston's code
% last modified LH & JFPJA & Bustardo 16/04/05

if (normalise)
    [J, C, R] = RAC_normaliser (J, mat_ref, rmask, 0); % apply RAC normilization to the matrix
end

% Inversion
[U, K, V] = svd(J*J'); %changed from nA because this program reads in a normalised sensitivity matrix

if isempty (no_svd)
    no_svd = length(find(diag(K)/K(1,1)>noise));
end
if PLOT
    hand=figure;
    ax_hand=axes('Parent',hand);
    semilogy(ax_hand,diag(K)/K(1,1));
    title ('SVs')
    xlabel('channels')
    ylabel('normalised SV')
end

num_meas = length(K);
% If singular value < tolerance, replace by 1e-6
k = spdiags(K);
% compensation for truncted energy 
%compensate = sum(k)/(sum(k(1:no_svd)));

newk = k;
newk(find(k<tolerance)) = 1e-6;
newk(no_svd+1:end) = 1e-6;

invk = 1./newk; % diagonal matrix inversion
invk(find(invk == 1e6)) = 0;
invK = spdiags(invk, 0, num_meas, num_meas);

% Construct inverted matrix Ainv
AA = V*invK*U';

if (normalise)
    ROI = find(rmask);
    if f_support

        [Sp,sp] = supporter_5('e', vtx, tri, []);
        %     if length(J)==length(vtx)
        %         [sp] = cond_shrinker2(tri, vtx, sp);
        %     end
        sp = sp/mean(sp);
        Sp = spdiags(sp(ROI), 0,length(ROI),length(ROI));
    else
        Sp = speye(length(ROI));
    end
    % Jinv = Sp*C*J'*AA; code is not memory efficient%
    %therefore do it in steps BP 05-10-11
    Jinv_t  =   Sp*C;
    clear Sp C
    J       =   J';
    Jinv_tt =   Jinv_t*J;
    clear Jinv_t J
    Jinv    =   Jinv_tt*AA;
    clear Jinv_tt AA
else
    Jinv = J'*AA;
end