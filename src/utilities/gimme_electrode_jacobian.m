function [elec_jac] = gimme_electrode_jacobian (Mesh, coords, u_drive, u_meas, areaORposition, z_c, df)
% GIMME_ELECTRODE_JACOBIAN returns jacobian based on electrode area or
% position.
%   MESH - should have Nodes_faces, Faces and pos defined on it
%   COORDS - matrix of the surface coordinate system
%   ORIGIN - array with x and y coordinate of origin of surface coordinate system
%   U_DRIVE - the potentials resulting from the unique driving currents with the last N entries being the electrode voltages
%   U_MEAS - the potentials resulting from the unique measurement currents with the last N entries being the electrode voltages
%   AREAORPOSITION - if 0 area-Jacobian, if 1 x-position-Jacobian, if 2 y-position-Jacobian
%   Z_C - contact impedance of the electrodes
%   SRF_AREA - surface areas of the electrodes

%% All we need are the surface informations of the Mesh, the protocol and corresponding surface and electrode potentials
vtx = Mesh.vtx; % all nodes of the mesh
srf = Mesh.srf; % all surface triangles
pos = Mesh.pos; % all electrode positions

[m,n]=size(coords);
coordinates = zeros(m*n/3,3);
for i = 1:n/3
    coordinates(1+(i-1)*m:i*m,:) = coords(:,1+(i-1)*3:i*3);
end
% coord0 = (origin(2)-1)/3*m + origin(1);


%==========================================================================%
% CHANGE POSITION OF ELECTRODE 30
%size_of_change = -70e-4 + 11.32e-4 + 12.26e-4 + 19.27e-4 + 9.72e-4 + 6.59e-4 + 8.07e-4 + 5.52e-4 + 8.45e-4 + 4.47e-4 + 4.16e-4 + 4.33e-4 + 5.81e-4 + 2.44e-4...
%                + 2.44e-4 + 3.22e-4 + 3.22e-4 + 3.11e-4 + 2.78e-4; %these are with filtered Jacobian jac(<1e-5) = inf
%size_of_change = size_of_change + 2.78e-3 - 0.47e-3; %from here on I used a step width of 1e-2 for faster convergence. ( i.e. (1./test(:,30))'*diff * -1e-3 )
%size_of_change = -70e-4 + 1.317e-3 + 1.631e-3 + 0.676e-3 + 1.236e-3 + 1.336e-3 + 0.521e-3 + 0.811e-3 + 0.409e-3; %stepsize is 1e-2 => inv(elec_jac(:,30)'*elec_jac(:,30)) * elec_jac(:,30)' * diff * -1e-2
%size_of_change = size_of_change + 7.722e-3 - 1.675e-3; %from here on I used a step width of 1e-1 for faster convergence.
%size_of_change = -70e-4 + 5.75e-4 + 45.83e-4 - 6.09e-4 - 5*4.15e-4 + 68.02e-4 - 3*3.07e-4 - 3*1.69e-4 - 5*4.98e-4; %this is for electrode 1 (very coarse). stepsize is 1e-2
%size_of_change = -70e-4 + 181.0e-4 - 56.09e-4 + 23.36e-4 - 15.02e-4 USW; %adjacent injection, step size is 1e-1
%size_of_change = 40e-4 - 17e-4 - 17e-4 - 8.5e-4 - 9.5e-4 - 5.3e-4 - 5.1e-4 - 5.1e-4 + -5e-4 - 5e-4; % electrode 31 which is very coarsely refined. Step size is 1e-2
%size_of_change = 0e-4;% - 15.5e-4 - 20.7e-4 - 7.3e-4 - 8.2e-4 - 4.9e-4 - 4.0e-4 - 4.0e-4 - 5*1.2e-4 - 5*0.9e-4;
%size_of_change = -70e-4 + 202.8e-4 - 107.4e-4 + 47.8e-4 - 2e-4; %electrode 1 on mesh with refined electrode 30. stepsize is 1e-2. protocol is normal 10 lines.
size_of_change = 0e-4;
electrode_nr = 31;

dist=(sum((coordinates-repmat(pos(electrode_nr,:),size(coordinates,1),1)).^2,2)).^0.5;
[~,idx]=min(dist);

if areaORposition == 1 % jacobian with respect to x-coordinate
    if dist(idx+m)<dist(idx-m) % electrode is closer to x+1 than x-1
        vector = coordinates(idx+m,:)-coordinates(idx,:);
        vector = vector/norm(vector)*size_of_change;
    else
        vector = coordinates(idx,:)-coordinates(idx-m,:);
        vector = vector/norm(vector)*size_of_change;
    end
else % jacobian with respect to y-coordinate
    if dist(idx+1)<dist(idx-1) % electrode is closer to y+1 than y-1
        vector = coordinates(idx+1,:)-coordinates(idx,:);
        vector = vector/norm(vector)*size_of_change;
    else
        vector = coordinates(idx,:)-coordinates(idx-1,:);
        vector = vector/norm(vector)*size_of_change;
    end
end

pos(electrode_nr,:) = pos(electrode_nr,:) + vector;
%=======================================================================%


elec_diam = 7e-3*ones(size(pos,1),1);
%elec_diam(30)  = 15.00e-3 - 2.86e-3 - 1.32e-3 - 0.67e-3 USW; %adjacent injection. step size is 10 to increase convergence speed.
%elec_diam(30) = 15.00e-3 - 2.899e-3 - 1.387e-3 - 0.991e-3 + 0.115e-3 - 0.951e-3 - 0.445e-3 + 0.910e-3 - 0.788e-3 - 0.048e-3; %step size 10, because changes are so small... Fourth iteration goes in wrong direction
%elec_diam(30) = elec_diam(30) - 0.628e-3; %increased step size to 100, to converge faster
[elec,~,srf_area,~] = set_electrodes_var_3(vtx(:,1:3),srf,pos,elec_diam,'s');
hold on;

% template Jacobian matrices
jacobians = {};

for l=1:size(pos,1)
    jac = sparse([],[],[],size(u_drive,1),size(u_drive,1),0);
    
    %% if Jacobian with respect to electrode position, compute the vector field
    if areaORposition>0
        dist=(sum((coordinates-repmat(pos(l,:),size(coordinates,1),1)).^2,2)).^0.5;
        [~,idx]=min(dist);
        %y_pos = mod(idx,m) - mod(coord0,m);
        %x_pos = (coord0-(idx-y_pos))/m;
        %disp(['electrode ', num2str(l), ' has coord ', num2str(coordinates(idx,:))]);
        %disp([num2str(pos(l,:)), ' ,  ', num2str(coordinates(idx-1,:)), ' ,  ' , num2str(coordinates(idx+1,:))]);
        if areaORposition == 1 % jacobian with respect to x-coordinate
            if dist(idx+m)<dist(idx-m) % electrode is closer to x+1 than x-1
                vector = coordinates(idx+m,:)-coordinates(idx,:);
                vector = vector/norm(vector);
            else
                vector = coordinates(idx,:)-coordinates(idx-m,:);
                vector = vector/norm(vector);
            end
        else % jacobian with respect to y-coordinate
            if dist(idx+1)<dist(idx-1) % electrode is closer to y+1 than y-1
                vector = coordinates(idx+1,:)-coordinates(idx,:);
                vector = vector/norm(vector);
            else
                vector = coordinates(idx,:)-coordinates(idx-1,:);
                vector = vector/norm(vector);
            end
        end
        quiver3(pos(l,1),pos(l,2),pos(l,3), vector(1)/100, vector(2)/100, vector(3)/100, 'r');
        hold on;
    end
    
    %% find the electrode boundary lines
    elec_srf = reshape(elec(l,:),3,[])';
    elec_srf = elec_srf(1:end-length(find(elec_srf==0))/3,:); % remove all zeros at the end of elec_srf
    S = [elec_srf(:,[1 2]); elec_srf(:,[1 3]); elec_srf(:,[2 3])];
    N = sort(S,2);
    M = sortrows(N);

    i = 1;
    inc = 1;
    bndry = [];

    while i <= size(M,1);
        if i==size(M,1)
            bndry(inc,:) = M(i,:);
            break;
        end
        ithrow = M(i,:);
        jthrow = M(i+1,:); 
        if ithrow == jthrow 
            i = i + 2;
        else
            % put in boundary node list
            bndry(inc,:) = M(i,:);
            inc = inc + 1;
            i = i + 1;
        end
    end
    
    %% Compute the Jacobian for each electrode boundary segment separately
    for j = 1:size(bndry,1)
        % identify the triangle this side belongs to
        facet = elec_srf(sum(ismember(elec_srf,bndry(j,:)),2)==2,:);
        
        % find outward normal vector for this side
        inner_node = vtx(facet(~ismember(facet,bndry(j,:))),:);
        outward_vector = (vtx(bndry(j,1),:)-inner_node) + (vtx(bndry(j,2),:)-inner_node);
        side_vector = vtx(bndry(j,1),:) - vtx(bndry(j,2),:);
        perpendicular = outward_vector - side_vector*dot(side_vector,outward_vector);
        normal = perpendicular/norm(perpendicular);
        
        % temporary: find opposite triangle to improve outward normal
        facets = srf(sum(ismember(srf,bndry(j,:)),2)==2,:);
        if size(facets,1)~=2
            disp('houston, we have a problem');
        else
            facet_neighbour = facets(sum(ismember(facets,facet),2)==2,:);
        end
        
        % vector field for Jacobian with respect to electrode radius
        if areaORposition == 0
            vector = (vtx(bndry(j,1),:)-pos(l,:)) + (vtx(bndry(j,2),:)-pos(l,:));
            vector = vector/norm(vector);
            quiver3(pos(l,1),pos(l,2),pos(l,3), vector(1)/100, vector(2)/100, vector(3)/100, 'r');
            hold on;
        end
        
        % setup sparse jacobian template for this boundary segment
        multiplier = - 1/(z_c*srf_area(l)) * dot(normal,vector) * norm(side_vector);
        jac( [bndry(j,1),bndry(j,2),size(u_drive,1)-size(pos,1)+l] , [bndry(j,1),bndry(j,2),size(u_drive,1)-size(pos,1)+l] ) = ...
            jac( [bndry(j,1),bndry(j,2),size(u_drive,1)-size(pos,1)+l] , [bndry(j,1),bndry(j,2),size(u_drive,1)-size(pos,1)+l] ) + ...
            multiplier* [  1/3 ,  1/6 , -1/2; ...
                          1/6 ,  1/3 , -1/2; ...
                         -1/2 , -1/2 ,  1 ]; % this matrix can for instance be gained using Gauss-Lobato quadrature
                                             % and the points (U_d-u_d_0),
                                             % (U_d-u_d_1), (U_m-u_m_0), (U_m-u_m_1) with midpoints
                                             
%         if (l==5)
%             %disp(['side coord 1: ', num2str(vtx(bndry(j,1),:)) , '\n']);
%             %disp(['side coord 2: ', num2str(vtx(bndry(j,2),:)) , '\n']);
%             %disp(['inner coord: ', num2str(inner_node) , '\n']);
%             disp(['outward vector: ', num2str(outward_vector) , '\n']);
%             disp(['side vector: ', num2str(side_vector) , '\n']);
%             disp(['perpendicular vector: ', num2str(perpendicular) , '\n']);
%             disp(['normal vector: ', num2str(normal) , '\n']);
%             disp(['field vector: ', num2str(vector) , '\n']);
%             disp(['multiplier: ', num2str(multiplier) , '\n']);
%         end
    end
    jacobians{l} = jac;
end

%% we now have the sparse template Jacobian matrices, which can be multiplied with the drive and measurement fields to get the corresponding Jacobian entries:
% elec_jac = zeros(size(u_drive,2)*size(u_meas,2),size(pos,1));
% for drive=1:size(u_drive,2) % iterate through all drive patterns
%     for meas=1:size(u_meas,2)
%         elec_jac((drive-1)*size(u_meas,2)+meas,:) = cellfun(@(a) dot(a*u_drive(:,drive),u_meas(:,meas)), jacobians); % multiply left and right for each electrode
%     end
% end

%% This uses the weird SuperSolver Format
elec_jac = zeros(sum(df),size(pos,1));
for drive = 1:size(u_drive,2)
    for meas = 1:df(drive)
        elec_jac(sum(df(1:drive-1))+meas,:) = cellfun(@(a) dot(a*u_drive(:,drive), u_meas(:,sum(df(1:drive-1)) + meas)), jacobians);
    end
end

end
