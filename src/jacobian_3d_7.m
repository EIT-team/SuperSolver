function [J] = jacobian_3d_7(V, v_f, Ela)
% Usage:[J] = jacobian_3d_7(V, v_f, Ela,);
%
% General:
% Calculates Jacobian (sensitivity) matrix at mat_ref
%
% ver. 6 allow construction of anisotropic Jacobian
%
% Input:
%  V - current field {n + no. of electrodes x no. of unique injection patterns}
%  Ela - normalised elements volumes {3*k x 3*k}
%  v_f - measurement field {n + no. of electrodes x no. of measurements}
%
% Output:
% J - Jacobian {no. of measurements x k} or for anisotropic case {no. of measurements x 6k} 
%------------------------------------------------------------------------------------------------------------------------
global Fem Fwd Mesh

if ~isempty(Mesh)
    [nv, dimen] = size(Mesh.vtx);  % number of vertices and dimension
    ns = length(Mesh.tri);  % number of elements
    Ela = sparse( (1:dimen*ns), (1:dimen*ns),kron( Mesh.support.',ones(1,dimen)) );
end

if sum(Fem.df) ~= size(v_f,2);
    error('Mismatched data input');
end

%Select the part referring to the interior nodes
V = Fwd.current_field(1:nv,:);
v_f = Fwd.measurement_field(1:nv,:);

% pre- allocate memory for jacobian construction
J = zeros(sum(Fem.df),ns);

cnt = 0;

% check whether to construct isotropic or anisotropic Jacobian
if Fem.isotropy
    for p = 1:size(V,2)
        DV =  Mesh.sgrad * V(:,p); %Gradient of the current fields
        for m = 1:Fem.df(p)
            Dvf = Mesh.sgrad * v_f(:,sum(Fem.df(1:p-1)) + m); %Gradient of the measurement fields
            Jrow_x3 = Dvf .* DV ;
            Jrow_u = Jrow_x3(1:3:end) + Jrow_x3(2:3:end) + Jrow_x3(3:3:end);
            cnt = cnt + 1;
            J(cnt,:) = -Jrow_u';
        end
    end
    % extract the support (elements volumes)
    Ela_d = spdiags(Ela(1:3:end,1:3:end));
    Ela_s = spdiags(Ela_d, 0, ns, ns);
    J = J * Ela_s; % multiply the integrand by the volume

else
    for p = 1:size(V,2)
        DV =  Mesh.sgrad * V(:,p); %Gradient of the current fields
        for m = 1:Fem.df(p)
            Dvf = Mesh.sgrad * v_f(:,sum(df(1:p-1)) + m); %Gradient of the measurement fields
            Jrow_x3 = Dvf .* DV ;
            cnt = cnt + 1;
            J(cnt,:) = -Jrow_x3';
        end
    end
    J = J * Ela; % multiply the integrand by the volume
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