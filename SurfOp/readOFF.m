function [pt,trg]=readOFF(fname,foutname)
%%%%%%%%  [pt trg]=ReadOFF(fname,foutname)
% readOFF.m
% This function is used to read the OFF file and return the point and
% triangle mesh.
%
% Input:
%   fname: the name of the OFF file
%   foutname: the name of the output file
%
% Output:
%   pt: the point mesh
%   trg: the triangle mesh
%
% Example:
%   [pt,trg]=readOFF('sphere.off','sphere_normal.off');

fid = fopen(fname,'r');
shape_type = fscanf(fid,'%s',1);

a = fscanf(fid,'%f',3);
num_pt  = a(1);
num_trg = a(2);

pt_coordinates = fscanf(fid,'%f',num_pt*3);
temp = fscanf(fid,'%f',num_trg*4);
pt = reshape(pt_coordinates,3,num_pt);
pt = pt';
temp = reshape(temp,4,num_trg);
temp = temp';
trg = temp(:,2:4) + 1;

fclose(fid);

if nargin == 2
    writeObjMesh2(pt,pt,trg,foutname);
    system(['./MeshNormal ' foutname]);
end
