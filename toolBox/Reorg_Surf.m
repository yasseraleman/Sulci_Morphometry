function [Surft] = Reorg_Surf(Surf);
%
% Syntax :
% [Surft] = Reorg_Surf(Surf);
%
% This scripts reorganize Surfaces in case of deleted faces. 
%
% Input Parameters:
%       Surf          :  Surface with deleted faces
%
% Output Parameters:
%       Surft         : Reorganized surface

% Related references:
%
%
% See also: 
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% April 9th 2012
% Version $1.0

Surft = Surf;
temp = unique(Surf.SurfData.faces) ;
neworder = [1:length(temp)]';
[a,t] = ismember(Surft.SurfData.faces,temp);
vert2delet = find(ismember([1:length(Surft.SurfData.vertices)],unique(Surft.SurfData.faces(:))) == 0);
Surft.SurfData.vertices(vert2delet,:) = [];
Surft.SurfData.faces = t;
if isfield(Surf.SurfData,'VertexNormals')
    Surft.SurfData.VertexNormals(vert2delet,:) = [];
end
if isfield(Surf.SurfData,'FaceVertexCData')
    Surft.SurfData.FaceVertexCData(vert2delet,:) = [];
end
if isfield(Surf,'Is')
    Surft.Is(vert2delet,:) = [];
end
if isfield(Surf,'Tri')
    [Tri] = Vert_Neib(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
    Temp = sum(Tri);
    Tri(:,Temp==0) = [];
    Surft.Tri = Tri;
end
return;