function [vtx_sorted,tri_sorted] = optimize6_c(vtx, tri, freal, PLOT)
% Usage:  [vtx_sorted,tri_sorted] = optimize6_c(vtx, tri, freal)
%
% General: 
% optimise connectivities mapping, by the following steps:
% 1. generates a sparse matrix of the nodal connections.
% 2. sorts the verteces ectors according to the sparse symmetric minimum degree ordering.
% 3. generates a sorted sparse matrix of the nodal connections.
% Actually the last phase is for the sake of visabilty only. those of you who rather spear wasting few 
% prescious seconds of their lifes, can omit this last section;
%
% ver 6. symmmd replaced by the more efficient symamd function. LH 14/4/5
%
% Input: 
% vtx - nodes vertices {n x 3}
% tri - elements as nodes indices {k x 4}
% freal - flag for complex numbers {scalar}
% PLOT - flag for plotting - 1 is on, 0 is off
%
% Output:
% vtx_sorted - sorted nodes vertices {n x 3}
% tri_sorted - sorted elements as nodes indices {k x 4}
%
% (!) impotant remark: be sure that the number of nodes is smaller than 2^16. or simply don't use the commands:
% uintri=uint16(tri); and uintri_sorted=uint16(tri_sorted); and then: uintri=tri; and uintri_sorted=tri_sorted;
%--------------------------------------------------------------------------------------------------------------------------------

% disp('optimising connectivity...');
% tic
nvtx = size(vtx,1);
ntri = size(tri,1);

if length(vtx)<2^16
uintri = uint16(tri); % since tri contains only indices we can save memory by using uint16, instead of double
% !!! impotant remark: be sure that the number of nodes is smaller than 2^16.
else
    uintri = tri;
end

combinations = nchoosek(1:4, 2); % contain all the different nodes connection combinations of an element 
combinations = [combinations;combinations(:,2) combinations(:,1)];

uintrix = uintri(:,combinations(:,1)'); % map all the x cooridantes of the sparse matrix newgr
uintriy = uintri(:,combinations(:,2)'); % map all the y cooridantes of the sparse matrix newgr

uintri_uniq = unique([uintrix(:) uintriy(:)],'rows'); % get rid of redandencies
dbtri_uniq=double(uintri_uniq);% since the sparse command can't work with uint16

gr = sparse(dbtri_uniq(:,1),dbtri_uniq(:,2),ones(size(uintri_uniq,1),1),nvtx,nvtx); %generates a sparse matrix according to the given nodal connections gr = gr + gr'; % complete the connectivity matrix with its symmetry
gr = gr + speye(nvtx);  %adds an identity sparse matrix since all nodes are self connected
% dur = toc; disp(toc)

clear uintri_uniq  uintrix uintriy% v_uintrix v_uintriy

if PLOT
    figure; 
    subplot(1,2,1);
    spy(gr);
    title('Connectivity before');
end
%--------------------------------------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------------------------------------
if freal
    nodesort = symamd(gr);
else
    nodesort = symrcm(gr);
end
clear gr

% disp('optimising nodes...');
tic
vtx_sorted = [];

i = 1:nvtx;
invsort(nodesort(1,i)) = i;
vtx_sorted(i,:) = vtx(nodesort(1,i),:);
% dur = toc; disp(toc)
%--------------------------------------------------------------------------------------------------------------------------------
% disp('optimising sorted connectivity, boundary vertices and param...');
tic
tri_sorted=[];

i = 1 : ntri;
tri_sorted(i,:) = invsort(tri(i,:));
% dur = toc; disp(toc)
%--------------------------------------------------------------------------------------------------------------------------------
% disp ('finally, generating new sorted connectivity matrix...')
tic

if length(vtx)<2^16
uintri_sorted = uint16(tri_sorted); % since tri_sorted contains only indices we can save memory by using uint16, instead of double
% !!! impotant remark: be sure that the number of nodes is smaller than 2^16.
else
    uintri_sorted = tri_sorted;
end

uintri_sortedx = uintri_sorted(:,combinations(:,1)'); % map all the x cooridantes of the sparse matrix newgr
uintri_sortedy = uintri_sorted(:,combinations(:,2)'); % map all the y cooridantes of the sparse matrix newgr

uintri_sorted_uniq = unique([uintri_sortedx(:) uintri_sortedy(:)],'rows'); % get rid of redandencies
dbtri_sorted_uniq=double(uintri_sorted_uniq);% since the sparse command can't work with uint16

newgr = sparse(dbtri_sorted_uniq(:,1),dbtri_sorted_uniq(:,2),ones(size(uintri_sorted_uniq,1),1),nvtx,nvtx); %generates a sparse matrix according to the given nodal connections
newgr = newgr + speye(nvtx);  %adds an identity sparse matrix since all nodes are self connected
% dur = toc; disp(toc)

if PLOT
    subplot(1,2,2)
    spy(newgr);title('Connectivitiy after');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This might become a part of the EIDORS suite
% Copyright (c) Lior Horesh 2004, EIT group, Medical Physics and Bioengineering, UCL, UK
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details
% EIDORS 3D version 2
% MATLAB Version 6.1.0.450 (R12.1) on PCWIN
% MATLAB License Number: 111672
% Operating System: Microsoft Windows Server 2003 Standard Edition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%