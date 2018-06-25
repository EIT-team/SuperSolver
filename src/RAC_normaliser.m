function [RAC, C, DR] = RAC_normaliser (J, mat_ref, rmask, PLOT)
% General: applys rows and columns element-based normalization to a matrix. Row normalization based on multiplication 
% of the sensitivity matrix by the conductivity vector, to produce V reference. Column normalization done 
% by taking the square root of the square norm.
% Written by LH on 22/9/03

% Input: J - Jacobian
%        mat_ref - conductivity
%        PLOT - flag for plotting, 0 - is off, 1 - is on
% Output: RAC - the normilized matrix.
%         R - sparse matrix of the row normalization.
%         C - sparse matrix of the column normalization.

global Sens
% C-Normalization
if ~isempty(find(ismember(Sens.norm,'c'))),
    c = 1./ (sqrt(mean(J(:,find(rmask)) .^2)))'; % or better be written as norm (J,'fro') The Frobenius-norm of matrix J, sqrt(sum(diag(J'*J)))
else
    c=ones(size(J,2),1);
end
C = spdiags(c,[0],length(c),length(c));

if PLOT
    figure, plot(abs(diag(C)));
    title ('Column normalisation') 
    xlabel('elements')
    ylabel('power')
end

% R-Normalization
if ~isempty(find(ismember(Sens.norm,'r'))),
    dr = J(:,find(rmask))*mat_ref(find(rmask)); % LH 150505
else
    dr=ones(size(J,1),1);
end
DR = spdiags(1./dr,0,length(dr),length(dr));
% DR = spdiags(1./dr,0,size(J,1),size(J,1)); LH 150505

if PLOT
    figure, plot(abs(diag(DR)));
    title ('Row normalisation') 
    xlabel('channels')
    ylabel('voltage')
end

% Applying normailzation
% RAC = DR*J*C; %RY code is not memory efficient
RA = DR*J(:,find(rmask));
clear J
RAC = RA*C;
clear RA