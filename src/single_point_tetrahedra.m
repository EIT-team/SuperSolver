function [tetrahedra_index]=single_point_tetrahedra(tri,vtx,point_of_interest,PLOT)

[centres]               =   find_tetrahedral_centres(tri,vtx);
dx                      =   centres(:,1)-point_of_interest(1);
dy                      =   centres(:,2)-point_of_interest(2);
dz                      =   centres(:,3)-point_of_interest(3);
vector_distance         =   sqrt((dx.^2+dy.^2+dz.^2));
[val tetrahedra_index]  =   min(vector_distance);

if PLOT
    figure;
    trimesh(dubs3_2(tri),vtx(:,1),vtx(:,2),vtx(:,3),'FaceColor','none','EdgeColor','b')
    daspect([1 1 1])
    hold on
    tetramesh(tri(tetrahedra_index,:),vtx,'FaceColor','r')
    
    
end