function varargout = Compute_Hull_from_Surface(Surf,varargin);
%
% Syntax :
% HullSurfMat = Compute_Hull_from_Surface(Surf,shFactor);
%
% This function computes hull surface from points coordinates.
%
% Input Parameters:
%   Surfa       : Surfaces files.
%   shFactor    : Shrinkage factor
%
% Output Parameters:
%   HullSurfMat       : Hull Surface.
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_oversurf Atlas_Surf
% Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 30th 2006
% Version $1.0

%% ====================== Checking input parameters ===================== %
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
    return;
else
    Surf = Surface_Checking(Surf);
    % Parameters
    opts.shfact = .95;  % No convex hull
    opts.outFilename = '';
end
if nargin == 2
    opts.shfact = varargin{1};
    if isnumeric(opts.shfact)
        error('The shrink factor must be a number between 0 and 1');
        return;
    end
    opts.outFilename = '';
end
if nargin == 3
    opts.shfact = varargin{1};
    opts.outFilename = varargin{2};
    if ~isnumeric(opts.shfact)
        error('The shrink factor must be a number between 0 and 1');
        return;
    end
end
%% =================== End of checking input parameters ================= %

%% ======================== Main Program ================================ %
HullSurfMat = Surf;

factors = linspace(opts.shfact,0,10);
Nf = length(factors);
for i = 1:Nf
    HullSurfMat = Surf;
    faces =   boundary(Surf.SurfData.vertices(:,1),Surf.SurfData.vertices(:,2),Surf.SurfData.vertices(:,3),factors(i));
    HullSurfMat.SurfData.faces = faces;
    HullSurfMat = Reorg_Surf(HullSurfMat);
    allVertex = HullSurfMat.SurfData.vertices;
    allFaces = HullSurfMat.SurfData.faces;
    allEdges = [allFaces(:,1) allFaces(:,2); allFaces(:,2) allFaces(:,3);allFaces(:,1) allFaces(:,3)] ; % Edges from the intersection faces
    allEdges = unique(sort(allEdges')','rows');
    Nedges = size(allEdges,1);
    Nvert = size(allVertex,1);
    Nfaces = size(allFaces,1);
    
    % Euler Characteristics
    eulerNumb = Nvert- Nedges + Nfaces;
    if eulerNumb == 2
        break;
    end
end
HullSurfMat = Surf;
HullSurfMat.SurfData.faces = faces;
HullSurfMat = Reorg_Surf(HullSurfMat);



HullSurfMat.SurfData = smoothpatch(HullSurfMat.SurfData,5);
HullSurfMat = Reorg_Surf(HullSurfMat);

if ~isempty(opts.outFilename);
    [OutFiles] = Save_Surf(HullSurfMat, opts.outFilename);
end
%========================End of main program==============================%
% Outputs;
varargout{1} = HullSurfMat;

return


