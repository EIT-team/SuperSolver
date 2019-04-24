function [bnd_v] = get_boundary_meas_2(Fwd,Fem,Mesh)
% Usage:  [bnd_v] = get_boundary_meas_2(V);
%
% General:
% Generates boundary voltages, for EIDORS based solutions.
%
% ver 2. assumes globals
% Input:
% V - approximated nodal potential distribution {n + no.electrodes x no. of unique current injections patterns}
% 
% Output:
% bnd_v - calculated boundary measurements {no. of measurements x 1}
%------------------------------------------------------------------------------------------------------------------------


V = Fwd.current_field;
    
curr_index = [];
bnd_v = zeros(1, size(Fem.prt,1));

for c = 1:length(Fem.df)
    curr_index = [curr_index; c * ones(Fem.df(c), 1)];
end

svtx = size(Mesh.vtx,1);
for p = 1:size(Fem.prt, 1)
    bnd_v(p) = (V((Fem.prt(p, 3) + svtx), curr_index(p)) - V((Fem.prt(p, 4) + svtx),curr_index(p)));
end

bnd_v = bnd_v';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This might become a part of the EIDORS suite
% Copyright (c) Lior Horesh & Rebecca Yerworth 2004, EIT group, Medical Physics and Bioengineering, UCL, UK
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details
% EIDORS 3D version 2
% MATLAB Version 6.1.0.450 (R12.1) on PCWIN
% MATLAB License Number: 111672
% Operating System: Microsoft Windows Server 2003 Standard Edition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%