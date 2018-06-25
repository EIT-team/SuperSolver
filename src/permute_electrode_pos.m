function [permuted_pos, size_of_change] = permute_electrode_pos(pos, std, coords)

    if (length(std)==1)
        size_of_change = std*randn(2*size(pos,1),1);
    else
        size_of_change = std;
    end
    permuted_pos = pos;

    [m,n]=size(coords);
    coordinates = zeros(m*n/3,3);
    for i = 1:n/3
        coordinates(1+(i-1)*m:i*m,:) = coords(:,1+(i-1)*3:i*3);
    end
    
    for electrode_nr = 1:size(pos,1)

        dist=(sum((coordinates-repmat(pos(electrode_nr,:),size(coordinates,1),1)).^2,2)).^0.5;
        [~,idx]=min(dist);

        % jacobian with respect to x-coordinate
        if dist(idx+m)<dist(idx-m) % electrode is closer to x+1 than x-1
            vector = coordinates(idx+m,:)-coordinates(idx,:);
            vector = vector/norm(vector)*size_of_change(2*electrode_nr-1);
        else
            vector = coordinates(idx,:)-coordinates(idx-m,:);
            vector = vector/norm(vector)*size_of_change(2*electrode_nr-1);
        end
        permuted_pos(electrode_nr,:) = permuted_pos(electrode_nr,:) + vector;
            
        % jacobian with respect to y-coordinate
        if dist(idx+1)<dist(idx-1) % electrode is closer to y+1 than y-1
            vector = coordinates(idx+1,:)-coordinates(idx,:);
            vector = vector/norm(vector)*size_of_change(2*electrode_nr);
        else
            vector = coordinates(idx,:)-coordinates(idx-1,:);
            vector = vector/norm(vector)*size_of_change(2*electrode_nr);
        end
        permuted_pos(electrode_nr,:) = permuted_pos(electrode_nr,:) + vector;
        
    end