function aroundpt=getaroundpt(pt,trg)
% getaroundpt.m
% This function is used to get the around points of a given point in a
% triangle mesh.
%
% Input:
%   pt: the point to get the around points of
%   trg: the triangle mesh
%
% Output:
%   aroundpt: the around points of the given point
if nargin<2
    disp('Parameter:  aroundpt=getaroundpt(pt,trg)');
    return
end

aroundpt=cell(length(pt),1);   

for i=1:length(trg)
    p1=trg(i,1);
    p2=trg(i,2);
    p3=trg(i,3);
    aroundpt{p1}=union(aroundpt{p1},[p1 p2 p3]);
    aroundpt{p2}=union(aroundpt{p2},[p1 p2 p3]);
    aroundpt{p3}=union(aroundpt{p3},[p1 p2 p3]);
end