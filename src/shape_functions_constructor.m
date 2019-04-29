function [D, sp] = shape_functions_constructor(vtx,tri)
% Usage: [D, sp] = shape_functions_constructor(vtx,tri);
%
% General: construct shape functions gradients and normalised volumes
%
%  Input:
%  vtx - vertices {n x 3}
%  tri - elements {k x 4}
%
% Output:
%  D - shape gradients of the shape functions over each element {3*k x n}
%  sp - elements volumes (support) {k x 1}
%
% 05/07/03 LH
%--------------------------------------------------------------------------
%if nargin == 0
%     vtx = Mesh.vtx;
%     tri = Mesh.tri;
% end

[nv, dimv] = size(vtx);  % number of vertices and dimension
ns = size(tri,1); % number of simplices

dimen2 = 2*dimv;

ils = 1:dimv*ns;
ilst(1:2:dimen2*ns) = ils;
ilst(2:2:dimen2*ns) = ils;

patt = 1:dimen2:ns*dimen2;

jlst(patt) = tri(:,1);
jlst(patt+1) = tri(:,2);
jlst(patt+2) = tri(:,1);
jlst(patt+3) = tri(:,3);
jlst(patt+4) = tri(:,1);
jlst(patt+5) = tri(:,4);

sls = ones(1,dimv*ns);
slst(1:2:dimen2*ns) = -sls;
slst(2:2:dimen2*ns) = sls;

% D0 is the matrix of the definitions of the gradients on elements
D0 = sparse(ilst,jlst,slst,dimv*ns, nv);

J1 = D0 * vtx(:,1);
J2 = D0 * vtx(:,2);
J3 = D0 * vtx(:,3);

r=1;
for w=1:3    
    JJ(r,w,:) = J1(w:dimv:ns*dimv);
    JJ(r+1,w,:) = J2(w:dimv:ns*dimv);
    JJ(r+2,w,:) = J3(w:dimv:ns*dimv);
end

V = squeeze(sum( [prod([JJ(1,1,:);JJ(2,2,:);JJ(3,3,:)]); prod([JJ(1,2,:);JJ(2,3,:);JJ(3,1,:)]);...
    prod([JJ(1,3,:);JJ(2,1,:);JJ(3,2,:)]); prod([-JJ(1,3,:);JJ(2,2,:);JJ(3,1,:)]);...
    prod([-JJ(1,1,:);JJ(2,3,:);JJ(3,2,:)]); prod([-JJ(1,2,:);JJ(2,1,:);JJ(3,3,:)])]));

ilst = kron((1:dimv*ns), ones(1,dimv));
jlst = zeros(1, ns*dimv^2);
for d = 1:dimv
    jlst(d:dimv:ns*dimv^2) = kron((d:dimv:dimv*ns),ones(1,dimv));
end
slst        =   zeros(1,ns*dimv^2);

pat         =   1:dimv^2:ns*dimv^2;

slst(pat)   =   squeeze(sum([prod([JJ(2,2,:) ; JJ(3,3,:)]); prod([-JJ(2,3,:) ; JJ(3,2,:)])])) ./ V;
slst(pat+1) =   squeeze(sum([prod([JJ(3,1,:) ; JJ(2,3,:)]); prod([-JJ(2,1,:) ; JJ(3,3,:)])])) ./ V;
slst(pat+2) =   squeeze(sum([prod([JJ(2,1,:) ; JJ(3,2,:)]); prod([-JJ(3,1,:) ; JJ(2,2,:)])])) ./ V;
slst(pat+3) =   squeeze(sum([prod([JJ(3,2,:) ; JJ(1,3,:)]); prod([-JJ(1,2,:) ; JJ(3,3,:)])])) ./ V;
slst(pat+4) =   squeeze(sum([prod([JJ(1,1,:) ; JJ(3,3,:)]); prod([-JJ(1,3,:) ; JJ(3,1,:)])])) ./ V;
slst(pat+5) =   squeeze(sum([prod([JJ(3,1,:) ; JJ(1,2,:)]); prod([-JJ(1,1,:) ; JJ(3,2,:)])])) ./ V;
slst(pat+6) =   squeeze(sum([prod([JJ(1,2,:) ; JJ(2,3,:)]); prod([-JJ(2,2,:) ; JJ(1,3,:)])])) ./ V;
slst(pat+7) =   squeeze(sum([prod([JJ(2,1,:) ; JJ(1,3,:)]); prod([-JJ(1,1,:) ; JJ(2,3,:)])])) ./ V;
slst(pat+8) =   squeeze(sum([prod([JJ(1,1,:) ; JJ(2,2,:)]); prod([-JJ(2,1,:) ; JJ(1,2,:)])])) ./ V;

LocJ        =   sparse(ilst,jlst,slst,dimv*ns,dimv*ns);
D           =   LocJ * D0;

sp          =   abs(V) / 6;

% elements support to be used later for the Jacobian matrix (does not include conductivity)
% Ela = sparse( (1:dimv*ns), (1:dimv*ns),kron( sp.',ones(1,dimv)) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) Lior Horesh 2004, EIT group, Medical Physics and Bioengineering, UCL, UK
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details
% MATLAB Version 6.1.0.450 (R12.1) on PCWIN
% MATLAB License Number: 111672
% Operating System: Microsoft Windows Server 2003 Standard Edition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%