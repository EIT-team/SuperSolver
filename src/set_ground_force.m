function [grd_ind] = set_ground_8(vtx,srf,grd_pos,cnts,PLOT)
% Usage: [grd_ind] = set_ground_8(vtx,srf,grd_pos,cnts)
%
% General:
% Pick ground index surface between the 2 central electrodes
%
% Input: 
% vtx - vertices matrix {n x 3}
% srf - external surfaces {no. of external surfaces x 3}
% grd_pos - position of the electrodes surfaces vertices {no. of electrodes x 3}
% cnts - centres of the surfaces {k x 3}
% PLOT - flag for plotting - 1 is on, 0 is off
%
% Output:
% grd_ind - index of the ground vertex {scalar}
%-------------------------------------------------------------------------------------------------------------

grd_ind = grd_pos;
grd_ind_v = grd_ind (ones(1,length(vtx)),:); % duplicate the electrode position, using tony's trick
dist_el = sum(((vtx-grd_ind_v).^2)');
[min_d, grd_ind] = min(dist_el);

if PLOT
    % plot on top of a surface mesh
    hold on;
    plot3(vtx(grd_ind,1),vtx(grd_ind,2),vtx(grd_ind,3),'.', ...
        'MarkerSize',32,'Color','b');

    text(vtx(grd_ind,1),vtx(grd_ind,2),vtx(grd_ind,3),'G','FontWeight','Bold','FontSize',16,'Color',[0.7 0.2 0.25]); % plots numbers on electrode pos
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