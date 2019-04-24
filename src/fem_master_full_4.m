function [E,Fem,Mesh] = fem_master_full_4(mat,Fem,Mesh)
% Usage: [E] = fem_master_full_4(mat);
%
% General:
% Builds up the system matrix based on the complete electrode model.
% accepts variant no. of surfaces per electrode
%
% Input:
% mat - conductivity vector {k x 1}
%
% Output:
% E - full rank system matrix based on the 3D complete electrode model {n + no. of electrodes x n + no. of electrodes}
%------------------------------------------------------------------------------------------------------------------------

[Ef,Fem,Mesh] = bld_master_full_3(mat,Fem,Mesh);
[E,Fem,Mesh] = ref_master_2(Ef,Fem,Mesh);

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