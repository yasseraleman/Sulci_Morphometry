function varargout = Area_Comp(varargin);
%
% Syntax :
% [At, Af] = Area_Comp(fv);
%
% This function computes the surface area for patch fv.
%
% Input Parameters:
%   fv               : Surface Patch.
%      fv.vertices   : Surface vertices
%      fv.faces      : Surface Triangles
%
% Output Parameters:
%   At       : Surfaces Area in cm^2.
%   Af       : Surfaces Area for each surface triangle in cm^2.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Atlas_Surf
% Plot_oversurf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% February 17th 2007
% Version $1.0

%% ===================== Checking Input Parameters ====================== %
if nargin == 0
    error('One Input is mandatory');
    return
elseif nargin == 1
    fv = varargin{1};
else
    error('Too many inputs');
    return
end
if isfield(fv,'SurfData')
    fv = fv.SurfData;
end
if ~isfield(fv,'vertices')|~isfield(fv,'faces')
    errordlg('Wrong Surface variable');
    return;
end
if nargout > 2
    errordlg('To Many Output Parameters');
    return;
end
%% ===================== End of Checking Input Parameters =============== %

%=========================Main program====================================%
d12 = sqrt(sum((fv.vertices(fv.faces(:,1),:) - fv.vertices(fv.faces(:,2),:)).^2,2));
d23 = sqrt(sum((fv.vertices(fv.faces(:,2),:) - fv.vertices(fv.faces(:,3),:)).^2,2));
d13 = sqrt(sum((fv.vertices(fv.faces(:,1),:) - fv.vertices(fv.faces(:,3),:)).^2,2));
per = (d12+d23+d13)/2;
Af = sqrt(per.*(per-d12).*(per-d23).*(per-d13))/100;% cm^2
At = sum(Af); % Area per each face
%========================End of main program==============================%

% Outputs
varargout{1} = At;
varargout{2} = Af;
return