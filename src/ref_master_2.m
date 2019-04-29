function [Er,Fem,Mesh] = ref_master_2(E,Fem,Mesh,sch)
% Usage:[Er] = ref_master_2(E,sch);
%
% General: applying reference to the system. Modifying the system matrix to preserve uniqueness
%
% ver 2. assumes globals
%
% Input:
% E - rank deficient by 1 system matrix
%
% Output:
% Er - full rank matrix
% sch - grounding scheme:
%        0 for grounding node gnd_ind
%        1 for grounding electrode gnd_ind
% gnd_ind - ground index
%
%==============================================================================%
% global Mesh Fem
[nv,jnk] = size(Mesh.vtx);

if nargin < 4
    sch = 0;
end

% ground a surface node
if sch == 0 
    Er = E;
    Er(Fem.gnd_ind,:)= 0;					%zeros(1,mas_c);
    Er(:,Fem.gnd_ind)= 0;					%zeros(mas_r,1);
    Er(Fem.gnd_ind,Fem.gnd_ind) = 1;

% ground one of the boundary electrodes
else
    Er = E;
    if Fem.gnd_ind > size(E,1) - nv
        error('Grounding electrode index out of range')
    end
    Er(nv+Fem.gnd_ind,:)= 0;					%zeros(1,mas_c);
    Er(:,nv+Fem.gnd_ind)= 0;					%zeros(mas_r,1);
    Er(nv+Fem.gnd_ind,nv+Fem.gnd_ind) = 1;
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