function [I, df] = set_currents_2(prt, pos, vtx, current)
% Usage: [I, df] = set_currents_2(prt, pos, vtx, current);
%
% General: 
% Construct currents matrix to be used as RHS of the forward solution
% ver. 2 enable controlling the uniform current injection magnitude
%
% Input:
%  prt - protocol (first 2 columns are injecting electrodes , last 2 columns are measurement electrodes) {no. of injection-measurement patterns x 4}
% pos - position of the electrodes surfaces vertices {no. of electrodes x 3}
% vtx - vertices {n x 1}
% current - current magintude in ampers {scalar}
%
% Output: 
% I - current vector for CEM {n + no. electrodes x no. of unique injection pairs}
% df - number of measurements combinations assigned to each current combination  {1 x no. of unique injection pairs}
%-------------------------------------------------------------------------------------------------------------

curr_prt = prt(:,1:2);
[curr_combn,i,j] = unique (curr_prt,'rows'); % taking only the unique combinations
unsorti = sort(i); % because we don't want matlab to sort our current combiations
curr_combn = curr_prt(unsorti,:); 

df = hist (j,1:length(i)); % counts how many times each current combination been used
df = df(j(unsorti)); % again to undo unique sorting

%Ib = zeros(size(pos,1),length(curr_combn));
% for i=1:length(curr_combn)
Ib = zeros(size(pos,1),size(curr_combn,1));
for i=1:size(curr_combn,1)
	Ib(curr_combn(i,1),i) = current;
	Ib(curr_combn(i,2),i) = -current;
end
I = [zeros(size(vtx,1),size(Ib,2)); Ib];

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