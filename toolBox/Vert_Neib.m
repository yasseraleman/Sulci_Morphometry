function [Tri] = Vert_Neib(faces,Nv,Nf);
%
% Syntax :
% [Tri] = Vert_Neib(faces,Nv,Nf);
%
% Computes the triangles that each point belongs. This function is just to
% call Vert_Neib.mex*
% 
% Input Parameters:
%   faces  : Mesh faces.
%   Nv     : Number of vertices.
%   Nf     : Number of faces.
%
% Output Parameters:
%   Tri          : Nv(Number of vertices)xN(maximun number of faces) matrix that contains the triangles
%                  that each point belongs.  
%
% Related references:
% 
%
% See also: Red_Surf Smooth_Surf Plot_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1st 2006
% Version $1.0

% c = accumarray(double(faces(:)),ones(size(faces(:),1),1));
% Tri = zeros(Nv,max(c)+2);
% Tri(:,1:2) =[[1:Nv]' c]; 
% for k = 1:Nv
%       ind = ismember(faces,k);
%       ind = find(sum(ind')>=1);
%       Tri(k,3:size(ind,2)+2) =ind;
% end
