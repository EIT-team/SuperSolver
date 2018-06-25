function [Ef] = bld_master_full_3(mat)
% Usage: [Ef] = bld_master_full_3(mat);
%
% General: system matrix assembling based on the complete electrode model 
% This function is called within fem_master_full_4
%
% ver .3 assumes globals
%
% Input:
% mat - conductivity vector {k x 1}
%
% vtx (global) - vertices {n x 3}
% tri (global) - simplices {k x 4}
% gnd_ind (global) - ground index (index of the vertex we consider as ground) {scalar}
% elec (global) - electrodes to elements assignment {no. electrodes x (max no. of faces representing electrode)*3}
% zc (global) - contact impedance {no. of electrodes x 1}
%
% Output:
% Ef - un-referenced system matrix {n x n}
%------------------------------------------------------------------------------------------------------------------------
global Mesh Fem Progress 

if nargin == 0 
    mat = Sol.current;
end
    
[vr, vc] = size(Mesh.vtx);
[sr, sc] = size(Mesh.tri);
[er, ec] = size(Fem.elec);

if length(mat) ~= sr
    error('Invalid conductivity information for this mesh');
end

if length(Fem.zc) == er
    %column vector zc with the contact impedances in [Ohms] is required

    [Ef] = bld_master_2(mat);

    vertos = sparse(zeros(vr,er));
    horos = sparse(zeros(er,vr+er));

    Ef = [Ef,vertos;horos];

    %Up to this point we have calculated the master matrix without the influence of contact impedance.
    while length(Fem.zc) ~= er
        disp(sprintf('Please enter the correct zc column vector with length: %d',er));
    end

    for q=1:er
        tang_area = 0;
        q_th_ele = Fem.elec(q,:);  % Select the row of nodes corresponding to the current electrode
        q_th_ele = nonzeros(q_th_ele)'; % LH 23/6/03 with advise of N. Polydorides, to allow variable no. of surface per electrode
        for w=1:3:length(q_th_ele)

            m = q_th_ele(w);
            n = q_th_ele(w+1);
            l = q_th_ele(w+2);

            % This way m & n nodes belong to the edge tangential to the electrode and also at the same simplex.
            % We will now evaluate the distance "tangential contact area" between m,n & l
            Are = triarea3d(Mesh.vtx([m n l],:));

            % area mnl
            cali_area = Are ./ (Fem.zc(q) * Fem.srf_area(q));  % coefficient for the area mnl
            tang_area = tang_area + cali_area;
            % Start modifying "expanding" the E master matrix

            Ef(m,vr+q) = Ef(m,vr+q) - cali_area/3 ; % Kv -> Ec  -> Vertical bar
            Ef(n,vr+q) = Ef(n,vr+q) - cali_area/3 ;
            Ef(l,vr+q) = Ef(l,vr+q) - cali_area/3 ;

            Ef(vr+q,m) = Ef(vr+q,m) - cali_area/3 ; % Kv' -> Ec' -> Horizontal bar
            Ef(vr+q,n) = Ef(vr+q,n) - cali_area/3 ;
            Ef(vr+q,l) = Ef(vr+q,l) - cali_area/3 ;

            Ef(m,m) = Ef(m,m) + cali_area/6; % Kz -> E -> Main bar
            Ef(m,n) = Ef(m,n) + cali_area/12;
            Ef(m,l) = Ef(m,l) + cali_area/12;

            Ef(n,m) = Ef(n,m) + cali_area/12;
            Ef(n,n) = Ef(n,n) + cali_area/6;
            Ef(n,l) = Ef(n,l) + cali_area/12;

            Ef(l,m) = Ef(l,m) + cali_area/12;
            Ef(l,n) = Ef(l,n) + cali_area/12;
            Ef(l,l) = Ef(l,l) + cali_area/6;

        end % dealing with this electrode

        Ef(vr+q,vr+q) = Ef(vr+q,vr+q) + tang_area;
    end %for the whole set of electrodes
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