function normal = vertnormal(mesh)
%% this is to approximate the out normal vector of each vertices by weighted average.
% we first compute the out unit normal of each faces of the triangles then,
% we take average for each vertices.

%% Initialize
node = mesh.node;
elem = mesh.elem;
N = size(node,1); % number of vertices
%NT = size(elem, 1); % number of triangles
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
% ve3 = node(elem(:,2),:)-node(elem(:,1),:);
cp = cross(ve1, ve2, 2); %
area =  sqrt(sum(cp.^2,2))/2; % area weighted average
%area = ones(NT,1); % simple average
normal=2*cp./area; % is there a bug?
nxArea = area.*normal(:,1);
nyArea = area.*normal(:,2);
nzArea = area.*normal(:,3);
patchArea = accumarray(elem(:),[area;area;area], [N 1]);
nxArea = accumarray(elem(:),[nxArea;nxArea;nxArea],[N 1]);
nyArea = accumarray(elem(:),[nyArea;nyArea;nyArea],[N 1]);
nzArea = accumarray(elem(:),[nzArea;nzArea;nzArea],[N 1]);
nx = nxArea./patchArea;
ny = nyArea./patchArea;
nz = nzArea./patchArea;
normal = [nx, ny, nz]./repmat(sqrt(nx.^2+ny.^2+nz.^2),1,3);