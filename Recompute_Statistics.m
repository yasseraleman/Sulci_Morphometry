function [varargout] = Recompute_Statistics(varargin);
%
% Syntax :
%  [IstatTable] = Recompute_Statistics(varargin);
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


%% ==================== End of Checking Input Parameters ================ %

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



%% ================ Computing Individual Nodes Metrics ================== %
outDirectory =  [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'SulcMorphometry'];
Varnames = [outDirectory filesep 'meshes' filesep opts.subjid '_SulcMorphometry_meshes_' opts.hemisphere '.mat'];
load(Varnames);

Nargs = size(argFiles,1);
for i = 1:Nargs
    argFile = deblank(argFiles(i,:));
    [metricsMatL, varNamesL] = Global_Sulci_Processing(argFile, sulcMetrics, reparmSulci, sulciLines, outDirectory, opts);
end
if Nargs == 1
    temp{1} = metricsMatL; clear metricsMatL
    metricsMatL = temp; clear temp
    temp{1} = varNamesL;  clear varNamesL
    varNamesL = temp; clear temp
    
end
clear sulcMetrics reparmSulci sulciLines;
%% ============= End of Computing Individual Nodes Metrics ============== %
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

%% ================ Computing Individual Nodes Metrics ================== %
outDirectory =  [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'SulcMorphometry'];
Varnames = [outDirectory filesep 'meshes' filesep opts.subjid '_SulcMorphometry_meshes_' opts.hemisphere '.mat'];
load(Varnames);

Nargs = size(argFiles,1);
for i = 1:Nargs
    argFile = deblank(argFiles(i,:));
    [metricsMatR, varNamesR] = Global_Sulci_Processing(argFile, sulcMetrics, reparmSulci, sulciLines, outDirectory, opts);
end
if Nargs == 1
    temp{1} = metricsMatR; clear metricsMatR
    metricsMatR = temp; clear temp
    temp{1} = varNamesR;  clear varNamesR
    varNamesR = temp; clear temp
    
end
clear sulcMetrics reparmSulci sulciLines;
%% ============= End of Computing Individual Nodes Metrics ============== %
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