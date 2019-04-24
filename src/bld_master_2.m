function [Ef,Mesh,Fem] = bld_master_2(mat,Mesh,Fem)
% Usage: [Ef] = bld_master_2(mat);
%
% General: Builds up the main compartment (GAP-SHUNT) of the system matrix for the complete electrode model. It is called within the function fem_master_full.
%
% Input:
% mat - reference admittivity at each element (in the complex case mat(i) = sigma(i) - epsilon(i)) {k x 1}
%
% ver 2. assumes globals and supports anisotropy
% 
% Output:
% Ef - unreferenced GAP based system matrix {n x n}
%
% 15/12/05 LH
%--------------------------------------------------------------------------
% global Mesh Fem

[n, dimv] = size(Mesh.vtx);  % number of vertices and dimension
k = size(Mesh.tri,1); % number of simplices

if nargin ==0 
    mat = Sol.current;
end

% generate sgradients of the shape functions over each element {3*k x n} and support (elements volumes) {k x 1}
[Mesh.sgrad, Mesh.support] = shape_functions_constructor(Mesh.vtx,Mesh.tri);

if length(mat) ~= length(Mesh.support)
    error('Some elements have not been assigned');
end;

%This is for the global conductance matrix (includes conductivity)
if Fem.isotropy
    Ela = sparse( (1:dimv*k), (1:dimv*k),kron( (mat.* Mesh.support).',ones(1,dimv)) );
else
    Ela1 = kron( Vols.',ones(1,dimv));
    Ela2 = [];
    for i=1:k
        Ela2 = sparse(blkdiag(Ela2,mat(3*i-2:3*i,3*i-2:3*i)));
    end
    Ela = diag(Ela1) .* Ela2;
end

Ef = Mesh.sgrad'*Ela*Mesh.sgrad;

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