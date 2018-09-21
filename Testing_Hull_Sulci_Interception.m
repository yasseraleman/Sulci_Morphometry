function Testing_Hull_Sulci_Interception


close all;
% % % % % % % % % %
opts.brainvisadir = '/media/COSAS/8-BrainVISADataBase-HCP';
opts.freesurferdir = '/media/HCPData/5-freesurfer_processing';
IdFile = '/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';
opts.bvsulcigroup = which('Brainvisa_nomenclature_sulci+STS.txt'); 
opts.hullsurf = 'freesurfer';
opts.boolparallel = 0;
if exist(deblank(IdFile),'file');
    Ids = char(textread(IdFile,'%s'));
else
    Ids = IdFile;
end
Ns = size(Ids, 1);
statFiles = '';
for i = 6:41
    opts.subjid = deblank(Ids(i,:));
      disp([' Processing Subject: ' opts.subjid  ' ====> ' num2str(i) ' of ' num2str(Ns)]);
        
        
        
%          sulcNames  = {'S.C._';'S.C.sylvian._';'F.C.M.ant._';'F.C.M.post._';'F.C.M.r.AMS.ant._';'F.P.O._';'S.T.s._';'F.Cal.ant.-Sc.Cal._'};
%         correctNodesRE = char(strcat(sulcNamesRE,repmat({'left'},[ length( sulcNamesRE) 1])));
        [~, ~, tempCorrectNodes] = Detecting_BrainVisa_Groups(which('Brainvisa_nomenclature_sulci+STS.txt'), 'left');
%         tempCorrectNodes = tempCorrectNodes(find(ismember(cellstr(tempCorrectNodes),correctNodesRE) == 0),:);
        for k = 1:size(tempCorrectNodes,1)
            ind = strfind(tempCorrectNodes(k,:),'_');
            sulcNames{k,1} = tempCorrectNodes(k,1:ind);
        end
        
%         sulcNames = {'S.C._'};
        
        outDirectory = [opts.brainvisadir  filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'SulcMorphometry'];
        
        
        %% =========================== Main Program ============================= %
%         Left Hemisphere
        Varnames = [opts.brainvisadir  filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'SulcMorphometry' filesep 'meshes' filesep opts.subjid '_SulcMorphometry_meshes_left.mat'];

        
        correctNodes = char(strcat(sulcNames,repmat({'left'},[ length( sulcNames) 1])));
        [sulchullIntercep] = Compute_Sulci_Hull_Interception(opts.freesurferdir, opts.brainvisadir, opts.subjid,'L', correctNodes,opts);
        load(Varnames);
        indsulci = find(cellfun(@isempty,reparmSulci));
        indhull  = find(cellfun(@isempty,sulchullIntercep));
    
        inddiffLeft = indhull(ismember(indhull,indsulci)==0);
        disp([' Left Hemisphere differences ' num2str(length(inddiffLeft))]);
         
        
        % Right Hemisphere
        Varnames = [opts.brainvisadir  filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'SulcMorphometry' filesep 'meshes' filesep opts.subjid '_SulcMorphometry_meshes_right.mat'];

        
        correctNodes = char(strcat(sulcNames,repmat({'right'},[ length( sulcNames) 1])));
        [sulchullIntercep] = Compute_Sulci_Hull_Interception(opts.freesurferdir, opts.brainvisadir, opts.subjid,'R', correctNodes,opts);
        load(Varnames);
        indsulci = find(cellfun(@isempty,reparmSulci));
        indhull  = find(cellfun(@isempty,sulchullIntercep));
        
        inddiffRight = indhull(ismember(indhull,indsulci)==0);
        disp([' Right Hemisphere differences ' num2str(length(inddiffRight))]);
        disp('  ' );
end
return;



function [sulchullIntercep] = Compute_Sulci_Hull_Interception(FreeSurferDatabaseDir, BrainVisaDatabaseDir, subjID, HemiChar, correctNodes,opts);


% ---- Creating Output Directories
opts.subjid = subjID;
%% ====================== General Processing ======================= %%
if ~exist([FreeSurferDatabaseDir filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii'],'file')
    cad = ['mri_convert -i ' [FreeSurferDatabaseDir filesep opts.subjid filesep 'mri' filesep 'ribbon.mgz'] ' -o ' [FreeSurferDatabaseDir filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii']];
    system(cad);
end
Vfs = spm_vol([FreeSurferDatabaseDir filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii']);

% Reading Talairach Transformation
Taltransf = [FreeSurferDatabaseDir filesep opts.subjid filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];

cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];



% Mandatory Files for the hemisphere
hullSFile   =    [FreeSurferDatabaseDir filesep opts.subjid filesep 'surf' filesep lower(HemiChar) 'h.wmhull'];
pialSFile   =    [FreeSurferDatabaseDir filesep opts.subjid filesep 'surf' filesep lower(HemiChar) 'h.pial'];
roiArgFile  =    [BrainVisaDatabaseDir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.arg' ];
switch HemiChar
    case 'L'
        hemiCad = 'left';
    case 'R'
        hemiCad = 'right';
end
Varnames = [BrainVisaDatabaseDir  filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'SulcMorphometry' filesep 'meshes' filesep opts.subjid '_SulcMorphometry_meshes_' hemiCad '.mat'];
load(Varnames);

%% %%%%%%%%%%%%%%%%%% For each Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
[ptharg, nmarg, ext]    = fileparts(roiArgFile);
SulcTmtktriFile = [ptharg filesep nmarg '.data' filesep 'aims_Tmtktri.mesh' ];
hullJunctPath   = [ptharg filesep nmarg '.data' filesep 'hull_junction.ima' ];
 
% Reading Hull
HullSurfMat = Read_Surface(hullSFile );
HullSurfMat.SurfData.vertices =HullSurfMat.SurfData.vertices+repmat(cras,[size(HullSurfMat.SurfData.vertices,1) 1]); % adding RAS center
HullSurfMat = Compute_Surface_Normals(HullSurfMat);

% Reading Hull
PialSurfMat = Read_Surface(pialSFile);
PialSurfMat.SurfData.vertices =PialSurfMat.SurfData.vertices+repmat(cras,[size(PialSurfMat.SurfData.vertices,1) 1]); % adding RAS center





% Reading Sulci
%% ========================== Reading Sulci Meshes ====================== %
% ---- Reading Tmtktri File
SulcTmtktri =     [BrainVisaDatabaseDir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.data' filesep 'aims_Tmtktri.gii' ];
SulcTmtktriFile = [BrainVisaDatabaseDir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar opts.subjid '_default_session_auto.data' filesep 'aims_Tmtktri.mesh' ];

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

Surfsulci = Read_Surface(SulcTmtktriFile);

% Converting to Brainvisa Coordinate system
Surfsulci = freesCS2brainvisaCS(Surfsulci,Vfs,'b2f'); % Converting to Brainvisa Coordinate system

% Read Arg Filoe
[NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(roiArgFile);
NodeNpoint = NodeNpoint(1:length(TmtktriIds));

    
%% ====================== Processing each sulcus NODE ================== %%
if isempty(correctNodes)
    correctNodes  = unique(SulcLabels,'rows');
end

failids = 0;
inds = find(ismember(cellstr(SulcLabels),cellstr(correctNodes)));
Nn = length(inds); % Number of Nodes in the sulci mesh
sulchullIntercep = cell(Nn,1);
varTemp = zeros(length(Surfsulci),1);
cont = 0;
for k = 1:Nn
    SulcSurfMat = Surfsulci(TmtktriIds(inds(k)));
    try
        [tempVar] = Intercept_Surface_with_Surface(HullSurfMat,  SulcSurfMat);
    catch
        tempVar = 'Nan';
    end
    if ~isempty(tempVar)
        varTemp(TmtktriIds(inds(k))) = 1;
    end
    if ~isempty(tempVar)&isempty(reparmSulci{TmtktriIds(inds(k))})
        try
            [Sulcmet, Surfo, SurfL] = Compute_Node_Metrics(SulcSurfMat, PialSurfMat, HullSurfMat);
        catch
            Sulcmet = '';
            Surfo = '';
            SurfL= '';
        end
        sulcMetrics{TmtktriIds(inds(k))} = Sulcmet;
        reparmSulci{TmtktriIds(inds(k))} = Surfo;
        sulciLines{TmtktriIds(inds(k))} = SurfL;
    end
    sulchullIntercep{TmtktriIds(inds(k))}  =  tempVar;
end
save(Varnames,'sulcMetrics','reparmSulci','sulciLines');
return;