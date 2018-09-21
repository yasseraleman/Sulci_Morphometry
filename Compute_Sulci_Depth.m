function varargout =  Compute_Sulci_Depth(varargin);

%
% Syntax :
%       [SulcVert] =  Compute_Sulci_Depth(Surf2Process,bottomlineCoords,toplineCoords,opts);
%
% This function computes depth lines that ride from top lines to the bottom
% lines of the sulcal median mesh.
%
%
% Input Parameters:
%      Surf2Process             : Sulcal Median Mesh (Matlab Format).
%      bottomlineCoords         : Coordinates of bottom lines.
%      toplineCoords            : Coordinates of top lines.
%      opts                     : Options:
%                                    - opts.ncurvp (Number of curve points).
%                                    - opts.sub_type (Kind of subvision applied. For more information type 'help perform_mesh_subdivision').
%                                    - opts.nsub (Mesh Subdivision iterations. For more information type 'help perform_mesh_subdivision').
%                                    - opts.verbose (Verbose boolean variable).
%
% Output Parameters:
%    
%      SulcVert                 : Depth profiles coordinates
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% February 13th 2015
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin<3 % the indispensable input arguments are not provided
    error('Three inputs are mandatory');
    return
else
    Surf2Process     = varargin{1};
    bottomlineCoords = varargin{2};
    toplineCoords    = varargin{3};
    
    Surf2Process = Surface_Checking(Surf2Process);
end
if nargin < 4
        opts.ncurvp = 10;             % Number of curve points.
        opts.sub_type = 'loop';       % Kind of subvision applied. For more information type 'help perform_mesh_subdivision'.
        opts.nsub = 3;                % Mesh Subdivision iterations. For more information type 'help perform_mesh_subdivision'.
        opts.verbose = 1;             % Verbose boolean variable.
else
    opts = varargin{4};
    if ~isfield(opts,'ncurvp')
        opts.ncurvp = 10;              % Number of curve points.
    end
    if ~isfield(opts,'sub_type')
        opts.sub_type = 'loop';       % Kind of subvision applied. For more information type 'help perform_mesh_subdivision'.
    end
    if ~isfield(opts,'nsub')
        opts.nsub = 3;                 % Mesh Subdivision iterations. For more information type 'help perform_mesh_subdivision'.
    end
    if ~isfield(opts,'verbose')
        opts.verbose = 1;              % Verbose boolean variable.
    end
end
if nargin > 5
    error('To many inputs');
    return;
end

if nargout > 2
    error('To many outputs');
    return;
end
%% ==================== End of Input parameters  ======================== %
if opts.verbose
    disp(' ');
    disp('Computing Sulcal Depth lines ... ');
    tic;
end

%% =============== Computing Top and Bottom Lines ======================= %

if opts.ncurvp > 10
    ndepth = round(opts.ncurvp/3); % Number of depth lines
else
    ndepth = 8;
end

opts.verb = 0;
[vertex,faces] = perform_mesh_subdivision(Surf2Process.SurfData.vertices',Surf2Process.SurfData.faces',opts.nsub,opts);
Surft.SurfData.faces = faces';
Surft.SurfData.vertices = vertex';

Graph = surface2graph(Surft);
SulcVert = [ 0 0 0 ];
for i = 1:size(bottomlineCoords,1)
    [~,locb] = min(sqrt((Surft.SurfData.vertices(:,1) - bottomlineCoords(i,1)).^2 + (Surft.SurfData.vertices(:,2) - bottomlineCoords(i,2)).^2 + (Surft.SurfData.vertices(:,3) - bottomlineCoords(i,3)).^2));
    
    [~,loct] = min(sqrt((Surft.SurfData.vertices(:,1) - toplineCoords(i,1)).^2 + (Surft.SurfData.vertices(:,2) - toplineCoords(i,2)).^2 + (Surft.SurfData.vertices(:,3) - toplineCoords(i,3)).^2));
    [D,PATH] = graphshortestpath(Graph, loct, locb);
    
    tempVar = [toplineCoords(i,:);(Surft.SurfData.vertices(PATH(1:end-1),:) + Surft.SurfData.vertices(PATH(2:end),:))/2;bottomlineCoords(i,:);];
    
    DepthCoords = fitCurveTo3DPts(tempVar, toplineCoords(i,:), bottomlineCoords(i,:),ndepth, 0); % Fitting a 3D curve
    
    if size(DepthCoords,1) ~= ndepth
        DepthCoords = Pts_fit_3Dline(tempVar, toplineCoords(i,:), bottomlineCoords(i,:),ndepth); % Fitting a 3D curve
    end
%     hold on;plot3(DepthCoords(:,1),DepthCoords(:,2),DepthCoords(:,3),'-c','Linewidth',2,'Markersize',20);hold on;plot3(DepthCoords(:,1),DepthCoords(:,2),DepthCoords(:,3),'.c','Linewidth',2,'Markersize',20);
    
    SulcVert = [SulcVert;DepthCoords];
end
SulcVert(1,:) = [];
for planes = 2:ndepth-1
    Temp = [planes:ndepth:(ndepth*(size(bottomlineCoords,1)-1)+planes)]';
    templine = SulcVert(Temp(1:end),:);
    
    templine = fitCurveTo3DPts(templine, templine(1,:), templine(end,:),size(bottomlineCoords,1), 0); % Fitting a 3D curve
    if size(templine,1) ~= ndepth
        templine = Pts_fit_3Dline(templine, templine(1,:), templine(end,:),size(bottomlineCoords,1)); % Fitting a 3D curve
    end
    templine = [smooth(templine(:,1)) smooth(templine(:,2)) smooth(templine(:,3))];

    SulcVert(Temp(1:end),:) = templine;
% %     hold on;plot3(templine(:,1),templine(:,2),templine(:,3),'-b','Linewidth',2,'Markersize',20);hold on;plot3(templine(:,1),templine(:,2),templine(:,3),'.c','Linewidth',2,'Markersize',20);
end

%% ============ End of Computing Top and Bottom Lines =========== %
% Outputs
varargout{1} = SulcVert;

if opts.verbose
    toc;
end
return;

% ======================= Internal Functions ============================ %
function [cPts] = Pts_fit_3Dline(XYZ, stPt, endPt, Npoints);
MAX_ITERATIONS = 5;
completedIterations = 0;

% Initialise control pts linearly between start/end anchors
cPts = interp1(0:1, [stPt; endPt], linspace(0,1,Npoints));

while completedIterations < MAX_ITERATIONS
    % Get the nearest-neighbour cntrl pt for each of the sample points
    sqDists = cellfun(@(x)sum(bsxfun(@minus, XYZ, x).^2,2), num2cell(cPts,2),'UniformOutput', false);
    [~, nnIdxs] = min(cat(2,sqDists{:}),[],2);
    % Keep the anchors, update inner cPts to the mean of their linked input pts
    for i = 2:Npoints-1
        cPts(i,:) = mean(XYZ(nnIdxs==i,:),1);
    end
    % Handle any cPts that didn't have linked pts so their mean became NaN
    goodIdxs = find(~isnan(cPts(:,1)));
    badIdxs = find(isnan(cPts(:,1)));
    cPts(badIdxs,:) = interp1(goodIdxs, cPts(goodIdxs,:), badIdxs);
    % Re-spread the control points out linearly
    cPtCumSumDists = cumsum([0; sqrt(sum(diff(cPts,1).^2,2))]);
    cPts = interp1(cPtCumSumDists, cPts, linspace(0, cPtCumSumDists(end), Npoints));
    completedIterations = completedIterations + 1;
end
return;