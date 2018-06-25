function [centres]=find_tetrahedral_centres(tri,vtx)

centres=(vtx(tri(:,1),:)+vtx(tri(:,2),:)+vtx(tri(:,3),:)+vtx(tri(:,4),:))/4;