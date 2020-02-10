function [elec, cnts, srf_area, center_of_mass] = set_electrodes_var_3(vtx, srf, pos, w, method, PlotFlag)
% Usage: [elec, cnts, srf_area] = set_electrodes_var_3(vtx, srf, pos, w, method, fig);
%
% General:
% Project electrodes on the mesh surface
%
% Input:
% vtx - vertices matrix {n x 3}
% srf - external surfaces {no. of external surfaces x 3}
% pos - position of the electrodes surfaces vertices {no. of electrodes x 3}
% w - electrode's diameter {scalar}
% method  - searching for elements below electrode surface: {char}
%           'n' - starting from the node which is closest to the electrode position
%           's' - starting from the element centre which is closest to the electrode position
%           'm' - allocating the surfaces whose most remote node is within radius w/2
%           'v' - allocating the surfaces whose centres are within radius w/2
%           from a mesh pointaccording to the centres
% fig - flag, if not exist surface mesh will be ploted {} (optional)
%
% Output:
% elec - boundary surfaces vertices indices, each row contains the vertices indices for one electrode,
% arranged in groups of 3's {no. electrodes x (max no. of faces representing electrode)*3}
% cnts - coordinates of the centre of each triangular boundary surface {no. of external surfaces x 3}
% srf_area - electrodes surface area {no. of electrodes x 1}
%-------------------------------------------------------------------------------------------------------------

if length(w) ==1
    w = w*ones(size(pos,1),1);
end

c_sels = cell (size(pos,1),1);
sels = [];
% vector containing the geometric center of each triangular  surface in x,y,z coordinates.
cnts = (vtx(srf(:,1),:) + vtx(srf(:,2),:) + vtx(srf(:,3),:))/3;

% center of mass vector
center_of_mass = zeros(size(pos,1),3);

if method == 's'
    for j = 1 : size(pos,1)
        
        % finds the closest surface to the electrode position
        el_pos = pos(j,:);
        el_pos_v = el_pos (ones(1,length(cnts)),:); % duplicate the electrode position, using tony's trick
        dist_el = sum(((cnts-el_pos_v).^2)');
        elec_srf = find(dist_el == min(dist_el)); % finds the index of the closest surface to the electrode
        
        %disp(['electrode ' num2str(j) ' is at ' num2str(cnts(elec_srf,:))]);
        %disp(num2str(cnts(elec_srf,:)))
        
        % finds the surfaces around the center of the electrode
        el_pos_srf = cnts(elec_srf,:);
        el_pos_srf_v = el_pos_srf (ones(1,length(cnts)),:); % duplicate the electrode position, using tony's trick
        dist_srf_el = sqrt(sum(((cnts-el_pos_srf_v).^2)'));
        
        % generates cell array, each entry contains the surfaces for a certain electrode
        c_sels{j} = find(dist_srf_el < ones(1,length(cnts))*w(j)/2); % finds the indices of the surfaces under the electrodes
        
        sels = [sels c_sels{j}]; % vector containing all the electrodes surfaces
        no_el_srf(j) = length(c_sels{j}); % stores the no of surfaces for each electrode
        
        srf_area (j) = 0; % initilize the total surface area of the electrode
        for jj = 1:no_el_srf(j)
            V = vtx(srf(c_sels{j}(jj),:)',:);
            [ta] = triarea3d(V); % calculate the the surface area of the surface
            srf_area(j) = srf_area(j) + ta;
            center_of_mass(j,:) = center_of_mass(j,:) + ta*cnts(c_sels{j}(jj),:);
        end
        
    end
    
elseif method == 'm'
    for j = 1 : size(pos,1)
        
        surfaces_vtx = [vtx(srf(:,1),:) , vtx(srf(:,2),:) , vtx(srf(:,3),:)];% surfaces vertices
        % finds the distance between the electrode location and the surfaces' vertices
        el_pos = pos(j,:);
        el_pos_v = repmat(el_pos,length(srf),3);% duplicate the electrode position
        dist_srf_el = max([sum((surfaces_vtx(:,1:3)-el_pos_v(:,1:3)).^2')' sum((surfaces_vtx(:,4:6)-el_pos_v(:,4:6)).^2')' sum((surfaces_vtx(:,7:9)-el_pos_v(:,7:9)).^2')']');
        
        % generates cell array, each entry contains the surfaces for a certain electrode
        c_sels{j} = find(dist_srf_el < ones(1,length(cnts))*w(j)/2); % finds the indices of the surfaces under the electrodes
        
        sels = [sels c_sels{j}]; % vector containing all the electrodes surfaces
        no_el_srf(j) = length(c_sels{j}); % stores the no of surfaces for each electrode
        
        srf_area (j) = 0; % initilize the total surface area of the electrode
        for jj = 1:no_el_srf(j)
            V = vtx(srf(sels(jj),:)',:);
            [ta] = triarea3d(V); % calculate the the surface area of the surface
            srf_area(j) = srf_area(j) + ta;
        end
        
    end
elseif method == 'v'
    for j = 1 : size(pos,1)
        
        
        el_pos = pos(j,:);
        % finds the surfaces around the center of the electrode
        el_pos_v = el_pos (ones(1,length(cnts)),:); % duplicate the electrode position, using tony's trick
        dist_srf_el = sqrt(sum(((cnts-el_pos_v).^2)'));
        
        % generates cell array, each entry contains the surfaces for a certain electrode
        c_sels{j} = find(dist_srf_el < ones(1,length(cnts))*w(j)/2); % finds the indices of the surfaces under with double radius area
        
        sels = [sels c_sels{j}]; % vector containing all the electrodes surfaces
        no_el_srf(j) = length(c_sels{j}); % stores the no of surfaces for each electrode
        
        srf_area (j) = 0; % initilize the total surface area of the electrode
        for jj = 1:no_el_srf(j)
            V = vtx(srf(sels(jj),:)',:);
            [ta] = triarea3d(V); % calculate the the surface area of the surface
            srf_area(j) = srf_area(j) + ta;
        end
        
    end
else
    for j = 1 : size(pos,1)
        
        % finds the closest surface to the electrode position
        el_pos = pos(j,:);
        el_pos_v = el_pos (ones(1,length(vtx)),:); % duplicate the electrode position, using tony's trick
        dist_el = sum(((vtx-el_pos_v).^2)');
        elec_srf = find(dist_el == min(dist_el)); % finds the index of the closest surface to the electrode
        
        % finds the surfaces around the center of the electrode
        el_pos_srf = vtx(elec_srf,:);
        el_pos_srf_v = el_pos_srf (ones(1,length(cnts)),:); % duplicate the electrode position, using tony's trick
        dist_srf_el = sqrt(sum(((cnts-el_pos_srf_v).^2)'));
        
        % generates cell array, each entry contains the surfaces for a certain electrode
        c_sels{j} = find(dist_srf_el < ones(1,length(cnts))*w(j)/2); % finds the indices of the surfaces under the electrodes
        
        sels = [sels c_sels{j}]; % vector containing all the electrodes surfaces
        no_el_srf(j) = length(c_sels{j}); % stores the no of surfaces for each electrode
        
        srf_area (j) = 0; % initilize the total surface area of the electrode
        for jj = 1:no_el_srf(j)
            V = vtx(srf(sels(jj),:)',:);
            [ta] = triarea3d(V); % calculate the the surface area of the surface
            srf_area(j) = srf_area(j) + ta;
        end
    end
end
%----------------------------------------------------------------------

elec = zeros(size(pos,1),max(no_el_srf)*3);
for j=1:size(pos,1)
    elec(j,1:no_el_srf(j)*3) = reshape(srf(c_sels{j},:)',1,no_el_srf(j)*3);
end
%----------------------------------------------------------------------
sels=sels';

center_of_mass = center_of_mass ./ (srf_area'*ones(1,3));

%-------------------------------plot the electrodes---------------------------------------
if PlotFlag
    figure, title('electrodes');
    trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3),'FaceAlpha',0.5);
    colormap([0 0 0]);
    daspect([1 1 1]);
    
    hold on;
    axis image;
    set(gcf,'Colormap',[0.6 0.7 0.9]);
    
    grid off
    view(3);
    hold on
    
    for u = 1:size(sels)
        paint_electrodes(sels(u),srf,vtx);
    end
    
    % add numbers to the electrodes
    el_caption = 1:size(pos,1);
    hold on;
    for i = 1:size(pos,1)
        text(pos(i,1),pos(i,2),pos(i,3),num2str(el_caption(i)),'FontWeight','Bold','FontSize',16,'Color',[0.4 0.2 0.65]); % plots numbers on electrode pos
    end
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