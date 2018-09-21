function [varargout] = Sulci_Nodal_Processing(varargin);
%
% Syntax :
%  [sulcMetrics, reparmSulci, sulciLines] = Sulci_Nodal_Processing(Surfsulci, PialSurfMat, HullSurfMat, opts);
%
% This script computes metrics for all the sulci nodes (BrainVisa Nodes).
%
% Input Parameters:
%       Surfsulci               : Sulci Surface in matlab format
%       PialSurfMat             : Hemispheric Pial Surface in Matlab Format
%       HullSurfMat             : Hemispheric Hull Surface in Matlab Format
%       opts                    : Options
%                              opts.maxsulcawidth (mm). Maximum Sulcal Width
%                              opts.mindist2hull (mm). Minimum distance to the hull surface mm
%                              opts.angthreshold (degrees). Minimum angle allowed between normals
%
% Output Parameters:
%       sulcMetrics             : Sulcus metrics
%       reparmSulci             : Reparametrized Sulcus
%        sulciLines             : Curves: Topline, Bottom line, Length,
%                                 Depth and width curves
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin < 4
    opts.maxsulcalwidth  = 12; % mm. Maximum Sulcal Width
    opts.maxsulcaldepth  = 70; % mm. Maximum Sulcal Depth
    opts.maxsulcallength = 200; % mm. Maximum Sulcal Length
    opts.boolparallel = 1;
else
    opts = varargin{4};
    if ~isfield(opts,'maxsulcalwidth')
        opts.maxsulcalwidth = 12;      % Maximum Sulcal Width
    end
    if ~isfield(opts,'maxsulcaldepth')
        opts.maxsulcaldepth = 70;           %  Maximum Sulcal Depth
    end
    if ~isfield(opts,'maxsulcallength')
        opts.maxsulcallength = 200;         % Maximum Sulcal Length
    end
    if ~isfield(opts,'ncurvp')
        opts.ncurvp = 40;     % Number of curve points
    end
    if ~isfield(opts,'boolparallel')
        opts.boolparallel = 1;     % Number of curve points
    end
end

if nargin < 3
    error('Three Inputs are needed');
    return
end
Surfsulci = varargin{1};
PialSurfMat = varargin{2};
HullSurfMat = varargin{3};
%% ====================== End of Input parameters  =======================%

% % % load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab.mat');



sulcMetrics = cell(length(Surfsulci),1);
reparmSulci = cell(length(Surfsulci),1);
sulciLines  = cell(length(Surfsulci),1);

Nsulc   = length(Surfsulci);
failids = 0;
tic;
% ========================= Parallel Running ============================ %
if opts.boolparallel
    clust = parcluster('local');
    numWorkers = clust.NumWorkers;
    try
        parpool(numWorkers-2);
    catch
        parpool(numWorkers-2);
    end
    parfor indsulc = 1:Nsulc
        SulcSurfMat = Surfsulci(indsulc);
        disp( ' ');
        disp( ' ');
        disp( ' ');
        disp([ ' Processing Sulcus ========>  ' num2str(indsulc) ' of ' num2str(Nsulc)]);

        try
            [Sulcmet, Surfo, SurfL] = Compute_Node_Metrics(SulcSurfMat, PialSurfMat, HullSurfMat);
        catch
            Sulcmet = '';
            Surfo = '';
            SurfL= '';
        end
        sulcMetrics{indsulc} =  Sulcmet;
        reparmSulci{indsulc} = Surfo;
        sulciLines{indsulc} = SurfL;
    end
    delete(gcp);
else % ====================== Secuential Running ======================== %
    for indsulc = 1:Nsulc
        SulcSurfMat = Surfsulci(indsulc);
        disp( ' ');
        disp( ' ');
        disp( ' ');
        disp([ ' Processing Sulcus ========>  ' num2str(indsulc) ' of ' num2str(Nsulc)]);
        try
            [Sulcmet, Surfo, SurfL] = Compute_Node_Metrics(SulcSurfMat, PialSurfMat, HullSurfMat);
        catch
            Sulcmet = '';
            Surfo = '';
            SurfL= '';
        end
        sulcMetrics{indsulc} =  Sulcmet;
        reparmSulci{indsulc} = Surfo;
        sulciLines{indsulc} = SurfL;
    end
end

% ========================= Correct and recompute ======================= %
for indsulc = 1:Nsulc
        SulcSurfMat = Surfsulci(indsulc);
    try
        [tempVar] = Intercept_Surface_with_Surface(HullSurfMat,  SulcSurfMat);
    catch
        tempVar = '';
    end
    if ~isempty(tempVar)&isempty(reparmSulci{indsulc})
        try
            [Sulcmet, Surfo, SurfL] = Compute_Node_Metrics(SulcSurfMat, PialSurfMat, HullSurfMat);
        catch
            Sulcmet = '';
            Surfo = '';
            SurfL= '';
        end
        sulcMetrics{indsulc} =  Sulcmet;
        reparmSulci{indsulc} = Surfo;
        sulciLines{indsulc} = SurfL;
    end
end
%% ================ End of Processing each sulcus NODE ================= %%
delete(gcp);

% Outputs
varargout{1} = sulcMetrics;
varargout{2} = reparmSulci;
varargout{3} = sulciLines;
toc;
%% ==================== End Processing each sulcus  ================ %%
return;
