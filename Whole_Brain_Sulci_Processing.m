function [varargout] = Whole_Brain_Sulci_Processing(varargin);
%
% Syntax :
%  [IstatTable] = Whole_Brain_Sulci_Processing(opts);
%
% This script computes metrics for all the sulci nodes in both hemispheres
% left and right.
%
% Input Parameters:
%       opts                    : Options:
%                              opts.maxsulcawidth % (mm). Maximum Sulcal Width
%                              opts.mindist2hull  % (mm). Minimum distance to the hull surface mm
%                              opts.angthreshold  % (degrees). Minimum angle allowed between normals
%                              opts.ncurvp        % Number of curve points
%                              opts.brainvisadir  % BrainVISA output
%                                                   directory
%                              opts.subjid        % Subject Id
%                              opts.hemisphere    % Hemisphere (left or right)
%                              opts.hullsurf      % Hull surface (freesurfer or brainvisa)
%                              opts.pialsurf      % Pial surface (freesurfer or brainvisa)
%                              opts.freesurferdir % FreeSurfer output
%                                                   directory
%
%
%
% Output Parameters:
%     IstatTable                : Statistics Output Filename
%
% See also: Global_Sulci_Processing Hemispheric_Sulci_Processing Sulci_Nodal_Processing  Compute_Node_Metrics
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0


%% ==================== Checking Input Parameters ======================= %
if nargin <1
    error('One input is mandatory');
    return;
end
opts = varargin{1};

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

if strcmp(opts.hullsurf,'freesurfer')
    if ~isfield(opts,'freesurferdir')     %(default BrainVisa _*Hemi.gii)
        error('Freesurferdir Directory is mandatory');
        return;
    else
        if ~exist(opts.freesurferdir, 'dir')
            error('Freesurferdir Directory does not exist');
            return;
        end
    end
    if ~exist([opts.freesurferdir filesep opts.subjid], 'dir')
        error('Subject Id does not exist inside the FreeSurfer Directory');
        return;
    end
end

%% ==================== End of Checking Input Parameters ================ %


%% ====================== General Processing ======================= %%
if strcmp(opts.hullsurf,'freesurfer')||strcmp(opts.pialsurf,'freesurfer')
    if ~exist([opts.freesurferdir filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii'],'file')
        mkdir([opts.freesurferdir filesep opts.subjid filesep 'tmp']);
        cad = ['mri_convert -i ' [opts.freesurferdir filesep opts.subjid filesep 'mri' filesep 'ribbon.mgz'] ' -o ' [opts.freesurferdir filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii']];
        system(cad);
    end
    Vfs = spm_vol([opts.freesurferdir filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii']);
    opts.freevol = Vfs; % FreeSurfer Volume
    delete([opts.freesurferdir filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii']);
end

if strcmp(opts.hullsurf,'freesurfer')
    if ~exist([opts.freesurferdir filesep opts.subjid filesep 'surf' filesep 'lh.wmhull'],'file')|~exist([opts.freesurferdir filesep opts.subjid filesep 'surf' filesep 'rh.wmhull'],'file'); % FreeSurfer Hull Volume)
%         [HullSurfaces] = Extract_WM_HullSurface(opts.freesurferdir, opts.subjid);
            inSurf  = [opts.freesurferdir filesep opts.subjid filesep 'surf' filesep 'lh.white'];
            outSurf = [opts.freesurferdir filesep opts.subjid filesep 'surf' filesep 'lh.wmhull'];
            HullSurfacesL = Compute_Hull_from_Surface(inSurf,1,outSurf);
            HullSurfaces = outSurf;
        
            inSurf  = [opts.freesurferdir filesep opts.subjid filesep 'surf' filesep 'rh.white'];
            outSurf = [opts.freesurferdir filesep opts.subjid filesep 'surf' filesep 'rh.wmhull'];
            HullSurfacesR = Compute_Hull_from_Surface(inSurf,1,outSurf);
            HullSurfaces = strvcat(HullSurfaces,outSurf);
    end
    
end
if strcmp(opts.hullsurf,'brainvisa')
    inSurf  = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep opts.subjid '_Lwhite.gii'];
    outSurf = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep opts.subjid '_Lwhite_hull.mesh'];
    HullSurfacesL = Compute_Hull_from_Surface(inSurf,1,outSurf);
    HullSurfaces = outSurf;

    
    inSurf  = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep opts.subjid '_Rwhite.gii'];
    outSurf = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep opts.subjid '_Rwhite_hull.mesh'];
    HullSurfacesR = Compute_Hull_from_Surface(inSurf,1,outSurf);
    HullSurfaces = strvcat(HullSurfaces,outSurf);

end
%% ====================== End of General Processing ================ %%

%% %%%%%%%%%%%%%%%%%% Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
opts.hemisphere = 'left';
HemiChar = 'L';
if exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ])&...
        exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto_lobar.arg' ])
    argFiles = strvcat([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ],...
        [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto_lobar.arg' ]);
elseif exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ])&...
        ~exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto_lobar.arg' ])
    argFiles = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ];
elseif ~exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ])&...
        exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto_lobar.arg' ])
    argFiles = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto_lobar.arg' ];
end


[metricsMatL, varNamesL] = Hemispheric_Sulci_Processing(argFiles, opts);


%% %%%%%%%%%%%%%%%%%% End of Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%% %%

%% %%%%%%%%%%%%%%%%%% Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
opts.hemisphere = 'right';
HemiChar = 'R';

if exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ])&...
        exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto_lobar.arg' ])
    argFiles = strvcat([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ],...
        [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto_lobar.arg' ]);
elseif exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ])&...
        ~exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto_lobar.arg' ])
    argFiles = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ];
elseif ~exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ])&...
        exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto_lobar.arg' ])
    argFiles = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto_lobar.arg' ];
end

[metricsMatR, varNamesR] = Hemispheric_Sulci_Processing(argFiles, opts);

%% %%%%%%%%%%%%%%%%%% End of Right Hemisphere %%%%%%%%%%%%%%%%%%%%%% %%


%% %%%%%%%%%%%%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
Nvar = 10;
for argi = 1:size(argFiles,1)
    [pth, nm, ext] = fileparts(deblank(argFiles(argi,:)));
    
    % ---- Detecting Save Cad
    if strcmp(nm(end-3:end),'auto')
        cad2Save = 'auto';
    elseif strcmp(nm(end-4:end),'lobar')
        cad2Save = 'lobar';
    else
        cad2Save = '';
    end
    
    switch opts.hullsurf
        case 'freesurfer'
            metricsDir = 'SulcMorphometry'; % Output Metrics Directory
        case 'brainvisa'
            metricsDir = 'SulcMorphometry_BVHull';
        otherwise
            if exist(opts.hullsurf,'file')
                metricsDir = 'CustomHulll';
            else
                error('Unreadable Hull Surface');
                return;
            end
    end
    
    % Reading
    L1 = cellstr(varNamesL{argi});
    R1 = cellstr(varNamesR{argi});
    Ls1 = reshape(L1,[Nvar length(L1)/Nvar]);
    Rs1 = reshape(R1,[Nvar length(R1)/Nvar]);
    TotNames = [Ls1;Rs1];
    TotNames = TotNames(:);
    
    
    L1 = metricsMatL{argi};
    R1 = metricsMatR{argi};
    Ls1 = reshape(L1,[Nvar length(L1)/Nvar]);
    Rs1 = reshape(R1,[Nvar length(R1)/Nvar]);
    TotVariab = [Ls1;Rs1];
    TotVariab = cellstr(num2str(TotVariab(:)));
    
    Temp = [TotNames';TotVariab'];
    cad2print = '';
    for vari = 1:size(Temp,2)
        cad2print = [cad2print char(Temp(:,vari)) [',';',']];
    end
    cad2print(:,end) = [];
    cad2print = [strvcat('Subjects',opts.subjid) [',';','] cad2print];
    warning off;
    mkdir([ opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep metricsDir filesep 'stats']);
    warning on;
    IstatTable = [ opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep metricsDir filesep 'stats' filesep opts.subjid '_sulcal_' cad2Save '_Stats.txt'];
    fidind = fopen(IstatTable,'wt');
    fprintf(fidind,'%s\n',cad2print(1,:));
    fprintf(fidind,'%s\n',cad2print(2,:));
    fclose(fidind);
    
end
varargout{1} =  IstatTable;
%% %%%%%%%%%%%%%%%%%%%%%%%% End of Statistics %%%%%%%%%%%%%%%%%%%%%% %%
return;