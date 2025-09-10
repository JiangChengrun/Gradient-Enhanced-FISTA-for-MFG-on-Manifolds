function viewVectF(Coord,Vect)
% This function is used to view the vector field
% Input: Coord = point coordinates
%        Vect = vector field
% Output: None


quiver3(Coord(:,1),Coord(:,2),Coord(:,3),Vect(:,1),Vect(:,2),Vect(:,3),...
    'Color','r','LineWidth',1,'AutoScale','on');
