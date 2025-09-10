function area = trianglearea(mesh)
%% COMPAREA Compute the area of surface triangular mesh
%   AREA = COMPAREA(MESH) return the area of the surface triangular meth.
%
%   Copyright (C) Hailong 
%   26/9/2022

%% Initalize
node = mesh.node;
elem = mesh.elem;

%% Compute the area
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
cp = cross(ve1, ve2, 2);
area =  sqrt(sum(cp.^2,2))/2;