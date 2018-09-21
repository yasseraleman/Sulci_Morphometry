function [varargout] = Hemispheric_Sulci_Processing(varargin);
%
% Syntax :
%  [metrics, metNames] = Hemispheric_Sulci_Processing(argFiles, opts, outDirectory);
%
% This script computes metrics for all the sulci nodes in the same hemisphere 
% and store group the results using different Arg Files (BrainVisa Nodes Organization).
%
% Input Parameters:
%       argFiles                : BrainVisa Arg Files for the same hemisphere.
%       opts                    : Options:
%                              opts.maxsulcawidth (mm). Maximum Sulcal Width
%                              opts.mindist2hull (mm). Minimum distance to the hull surface mm
%                              opts.angthreshold (degrees). Minimum angle allowed between normals
%                              opts.ncurvp = 40; % Number of curve points
%      outDirectory             : Output Directory.
%
%
% Output Parameters:
%     metrics                   : Cellarray where each cell containings the
%                                 obtained metrics for each Arg file.
%     metNames                  : Metrics Names.
%
% See also: Global_Sulci_Processing Sulci_Nodal_Processing  Compute_Node_Metrics
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0


%% ==================== Checking Input Parameters ======================= %
if nargin <2
    error('At least two inputs are mandatory');
    return;
end
argFiles = varargin{1};
opts = varargin{2};

 
if ~isfield(opts,'brainvisadir')     %(default BrainVisa _*Hemi.gii)
    error('BrainVisa Directory is mandatory');
    return;
else
    if ~exist(opts.brainvisadir, 'dir')
        error('BrainVisa Directory does not exist');
        return;
    end
end
if ~isfield(opts,'subjid')
    error('Subject ID is mandatory');
    return;
else
    if ~exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid], 'dir')
        error('Subject Id does not exist inside the BrainVisa Directory');
        return;
    end
end

if ~isfield(opts,'pialsurf')
    opts.pialsurf = 'brainvisa'; % Pial Surface file (default BrainVisa _*Hemi.gii)
end
if ~isfield(opts,'hullsurf')     %(default BrainVisa _*Hemi.gii)
    opts.hullsurf = 'brainvisa';
end
if ~isfield(opts,'pialsurf')
    opts.pialsurf = 'brainvisa'; % Pial Surface file (default BrainVisa _*Hemi.gii)
end
if ~isfield(opts,'hemisphere') % Hemisphere
    opts.hemisphere = 'left';
    warning('Left Hemisphere has been selected as default');
end
if ~isfield(opts,'maxsulcawidth') % mm. Maximum Sulcal Width
    opts.maxsulcawidth = 12;
end
if ~isfield(opts,'maxsulcaldepth')
    opts.maxsulcaldepth = 70;           %  Maximum Sulcal Depth
end
if ~isfield(opts,'maxsulcallength')
    opts.maxsulcallength = 200;         % Maximum Sulcal Length
end
if ~isfield(opts,'ncurvp') % Number of curve points
    opts.ncurvp = 40;
end

switch opts.hemisphere
    case 'left'
        HemiChar = 'L';
    case 'right'
        HemiChar = 'R';
end
%% ==================== End of Checking Input Parameters ================ %

%% ========================== Reading Hull Surface ====================== %

changesurf = 0;
if strcmp(opts.hullsurf,'freesurfer')|strcmp(opts.pialsurf,'freesurfer')
    changesurf = 1; % Boolean variable to convert surfaces from brainvisa space to freesurfer space
    if ~exist([opts.freesurferdir filesep filesep opts.subjid], 'dir')
        error('Subject Id does not exist inside the FreeSurfer Directory');
        return;
    end
    if isempty(opts.freevol)
        cad = ['mri_convert -i ' opts.freesurferdir  filesep opts.subjid filesep 'mri' filesep 'ribbon.mgz -o ' opts.freesurferdir  filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii'] ;
        system(cad);
        Vfs = spm_vol([opts.freesurferdir  filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii']);
        opts.freevol = Vfs;
        delete(Vfs.fname);
    end
end

switch opts.hullsurf
    case 'freesurfer'
        metricsDir = 'SulcMorphometry'; % Output Metrics Directory
        HullSurfMat = Read_Surface( [opts.freesurferdir filesep opts.subjid filesep 'surf' filesep lower(HemiChar) 'h.wmhull']);
        Taltransf = [opts.freesurferdir filesep opts.subjid filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
        
        % Reading Talairach Transformation
        cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
        HullSurfMat.SurfData.vertices =HullSurfMat.SurfData.vertices+repmat(cras,[size(HullSurfMat.SurfData.vertices,1) 1]); % adding RAS center
        %HullSurfMat = freesCS2brainvisaCS(HullSurfMat,opts.freevol,'f2b'); % Converting to FreeSurfer Coordinate system
    case 'brainvisa'
        metricsDir = 'SulcMorphometry_BVWHull';
        BVHullSurfMesh = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep opts.subjid '_' HemiChar 'white_hull.mesh'];
        HullSurfMat = Read_Surface(BVHullSurfMesh);
        if changesurf
            HullSurfMat = freesCS2brainvisaCS(HullSurfMat,opts.freevol,'b2f'); % Converting to FreeSurfer Coordinate system
        end
    case 'brainvisa_hemi'
        metricsDir = 'SulcMorphometry_BVHull';
        % Selecting BrainVisaHull
        BVHullSurf =     [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep opts.subjid '_' HemiChar 'hemi_hull.gii'];
        HullSurfMat = Read_Surface(BVHullSurf);
        if changesurf
            HullSurfMat = freesCS2brainvisaCS(HullSurfMat,opts.freevol,'b2f'); % Converting to FreeSurfer Coordinate system
        end
    otherwise
        if exist(opts.hullsurf,'file')
            HullSurfMat = Read_Surface(opts.hullsurf);
        else
            error('Unreadable Hull Surface');
            return;
        end
end

% Computing Hull Normals
HullSurfMat = Compute_Surface_Normals(HullSurfMat);

% Output directory
if nargin < 3
    outDirectory = [ opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep metricsDir];
else
    outDirectory = varargin{3};
end
%- ==================== End of Reading Hull Surface ===================== %

%% ===================== Reading  Pial Surface ========================== %
switch opts.pialsurf
    case 'brainvisa'
        % 1. ------- From Brainvisa --------- %
        HemiGII  = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep opts.subjid '_' upper(HemiChar) 'hemi.gii'];
% % %         HemiMesh = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep opts.subjid '_' upper(HemiChar) 'hemi.mesh'];
% % %         if ~exist(HemiMesh,'file')
% % %             cad = ['AimsFileConvert -i ' HemiGII ' -f GIS -o ' HemiMesh ' --verbose 0'];
% % %             system(cad);
            try
                PialSurfMat = Read_Surface(HemiGII);
            catch
                stempSurf = gifti(HemiGII);
                PialSurfMat.SurfData.vertices = stempSurf.vertices;
                PialSurfMat.SurfData.faces = stempSurf.faces;
            end
% % %         else
% % %             PialSurfMat = Read_Surface(HemiMesh);
% % %         end
        if changesurf
            PialSurfMat = freesCS2brainvisaCS(PialSurfMat,opts.freevol,'b2f'); % Converting to Brainvisa Coordinate system
        end
        % 2. ------- From Freesurfer --------- %
        %%%%%%%%%%%%%%%load('/media/COSAS/scripts/Sulcal_Processing/matlab.mat');
    case 'freesurfer'
        lhsurf = [opts.freesurferdir filesep opts.subjid filesep 'surf' filesep lower(HemiChar) 'h.pial'];
        % Reading freesurfer white surface
        PialSurfMat= Read_Surface(lhsurf);
        PialSurfMat.Name = 'LH.WHITE';
        Taltransf = [opts.freesurferdir filesep opts.subjid filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
        
        % Reading Talairach Transformation
        cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
        PialSurfMat.SurfData.vertices =PialSurfMat.SurfData.vertices+repmat(cras,[size(PialSurfMat.SurfData.vertices,1) 1]); % adding RAS center
end

PialSurfMat = Compute_Surface_Normals(PialSurfMat);

%% ==================== End of Reading Pial Surface ===================== %


%% ========================== Reading Sulci Meshes ====================== %
% ---- Reading Tmtktri File
SulcTmtktri =     [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.data' filesep 'aims_Tmtktri.gii' ];
SulcTmtktriFile = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.data' filesep 'aims_Tmtktri.mesh' ];

if ~exist(SulcTmtktriFile,'file')
    cad = ['AimsFileConvert -i ' SulcTmtktri ' -f GIS -o ' SulcTmtktriFile ' --verbose 0'];
    system(cad);
end

% Loading Sulci Mesh Surface
 Surf = Read_Surface(SulcTmtktriFile);
% Surf = Read_Surface(SulcTmtktri);


contdelsurf = 0;
scont = 0;
indel = 0;
for j = 1:length(Surf)
    if isempty(Surf(j).Tri)
        contdelsurf = contdelsurf + 1;
        indel(contdelsurf) = j;
    else
        scont = scont + 1;
        sizes(scont) = size(Surf(j).Tri,1);
    end
end
if sum(indel) > 0
    % Removing empty Surfaces
    Surf(indel) = [];
end

% Saving the new mesh surface
save_mesh(Surf, SulcTmtktriFile);
Surfsulci = Surf; clear Surf;
if changesurf
    Surfsulci = freesCS2brainvisaCS(Surfsulci,opts.freevol,'b2f'); % Converting to FreeSurfer Coordinate System
end
%% =================== End of Reading Sulci Meshes ====================== %

%% ================= Computing Individual Nodes Metrics ================= %
[sulcMetrics, reparmSulci, sulciLines] = Sulci_Nodal_Processing(Surfsulci, PialSurfMat, HullSurfMat, opts);
% load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/resmatlab.mat');
warning off;
mkdir([outDirectory filesep  'meshes']);
OutdirMesh = [outDirectory filesep 'meshes'];
%system(['rm -r ' outDirectory filesep 'meshes' filesep '*']);
warning on;
Varnames = [OutdirMesh filesep opts.subjid '_' metricsDir '_meshes_' opts.hemisphere '.mat'];
save(Varnames,'sulcMetrics','reparmSulci','sulciLines'); 
%load(Varnames);
%% ============= End of Computing Individual Nodes Metrics ============== %
Nargs = size(argFiles,1);
for i = 1:Nargs
    argFile = deblank(argFiles(i,:));
    [metricsMat, measNames] = Global_Sulci_Processing(argFile, sulcMetrics, reparmSulci, sulciLines, outDirectory, opts);
    metrics{i}  =  metricsMat;
    metNames{i} =  measNames;
end
clear sulcMetrics reparmSulci sulciLines;
varargout{1} = metrics;
varargout{2} = metNames;
return;



