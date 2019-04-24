function [j,cnts,sp] = compute_current_densities(mesh,v)
% Computes the current density on a mesh based on the potential
% distribution
%
% input: mesh       mesh with vtx and tri and mat_ref (struct)
%        v          potentials on all nodes of the mesh (vector)
%
% output: j         current density vector for each element
%         cnts      center of each element
%
% Markus Jehl, 13/02/2014
 
%%%% firstly, we have to compute the potential gradient on each element %%%%
[D, sp] = shape_functions_constructor(mesh.vtx,mesh.tri);
Dv = D*v;

% cnts = zeros(size(mesh.tri,1),3);
% Dv = zeros(size(mesh.tri,1),3);
% for i=1:size(mesh.tri,1)
%    [x_min,argmin] = min(mesh.vtx(mesh.tri(i,:),:));
%    [x_max,argmax] = max(mesh.vtx(mesh.tri(i,:),:));
%    v_min= v(mesh.tri(i,argmin))';
%    v_max= v(mesh.tri(i,argmax))';
%    Dv(i,:)=(v_max-v_min)./(x_max-x_min);
%    
%    cnts(i,:) = mean(mesh.vtx(mesh.tri(i,:),:), 1);
% end
% sigma = mesh.mat_ref;
% j = [Dv(:,1).*sigma, Dv(:,2).*sigma, Dv(:,3).*sigma];

%%%% then we multiply the gradients with the respective conductivity %%%%
sigma(1:3:3*size(mesh.mat_ref,1)) = mesh.mat_ref;
sigma(2:3:3*size(mesh.mat_ref,1)) = mesh.mat_ref;
sigma(3:3:3*size(mesh.mat_ref,1)) = mesh.mat_ref;
j = Dv.*sigma';
j = [j(1:3:end),j(2:3:end),j(3:3:end)];

cnts = zeros(size(mesh.tri,1),3);
for i=1:size(mesh.tri,1)
    cnts(i,:) = mean(mesh.vtx(mesh.tri(i,:),:),1);
end

% for i=1:3
%     [x_min,argmin] = min(reshape(mesh.vtx(mesh.tri(:),i),size(mesh.tri,1),4)');
%     [x_max,argmax] = max(reshape(mesh.vtx(mesh.tri(:),i),size(mesh.tri,1),4)');
%     v_min = v(mesh.tri(
% end