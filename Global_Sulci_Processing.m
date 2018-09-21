function [varargout] = Global_Sulci_Processing(varargin);
%
% Syntax :
%  [metricsMat, measNames] = Global_Sulci_Processing(argFile, sulcMetrics, reparmSulci, sulciLines, outDirectory, opts);
%
% This script computes metrics for all the sulci nodes (BrainVisa Nodes).
%
% Input Parameters:
%       argFile                 : BrainVisa Arg File.
%       sulcMetrics             : Sulcus metrics.
%       reparmSulci             : Reparametrized Sulcus.
%        sulciLines             : Curves: Topline, Bottom line, Length,
%                                 Depth and width curves.
%      outDirectory             : Output Directory.
%       opts                    : Options:
%                              opts.maxsulcawidth (mm). Maximum Sulcal Width
%                              opts.mindist2hull (mm). Minimum distance to the hull surface mm
%                              opts.angthreshold (degrees). Minimum angle allowed between normals
%                              opts.ncurvp = 40; % Number of curve points
%
%
% Output Parameters:
%
%
%
% See also: Sulci_Nodal_Processing  Compute_Node_Metrics
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin < 6
    opts.bvsulcigroup = which('Brainvisa_nomenclature_sulci+STS.txt');
else
    opts = varargin{6};
end
if ~isfield(opts,'bvsulcigroup') % Number of curve points
    opts.bvsulcigroup = which('Brainvisa_nomenclature_sulci+STS.txt');
end

if nargin < 5
    error('Five Inputs are needed');
    return
end
argFile = varargin{1};
sulcMetrics = varargin{2};
reparmSulci = varargin{3};
sulciLines = varargin{4};
outDirectory = varargin{5};
%% ====================== End of Input parameters  =========================%

% load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab_test_Results.mat')
% Loading argFile
[NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(argFile);
NodeNpoint = NodeNpoint(1:length(TmtktriIds));



%% ===================== Detecting BrainVisa Output ===================== %
[pth, nm, ext] = fileparts(argFile);

% ---- Detecting Hemisphere
HemiChar = nm(1);
ind = strfind(argFile,'subjects');

% ---- Detecting BrainVisa Database Directory
BrainVisaDatabaseDir = argFile(1:ind-2);


% ---- Detecting Subject ID
ind = strfind(nm,'_default_session_auto');
subjId = nm(2:ind-1);

% ---- Detecting Save Cad

if strcmp(nm(end-3:end),'auto')
    cad2Save = 'auto';
elseif strcmp(nm(end-4:end),'lobar')
    cad2Save = 'lobar';
else
    cad2Save = '';
end


% ---- Creating Output Directories
%     mkdir([outDirectory filesep  cad2Save '_meshes']);
%     OutdirMesh = [outDirectory filesep cad2Save '_meshes'];
%     system(['rm -r ' outDirectory filesep cad2Save '_meshes' filesep '*']);
warning off;
mkdir([outDirectory filesep cad2Save '_sulcalspams']);
OutdirSPAM = [outDirectory filesep cad2Save '_sulcalspams'];
% system(['rm -r ' outDirectory filesep cad2Save '_sulcalspams' filesep '*']);

%
mkdir([outDirectory filesep cad2Save '_sulcallength']);
OutdirLENGTH = [outDirectory filesep cad2Save '_sulcallength'];
% system(['rm -r ' outDirectory filesep cad2Save '_sulcallength' filesep '*']);
%
mkdir([outDirectory filesep cad2Save '_sulcaldepth']);
OutdirDEPTH = [outDirectory filesep cad2Save '_sulcaldepth'];
% system(['rm -r ' outDirectory filesep cad2Save '_sulcaldepth' filesep '*']);
warning on;
%% =============== End of  Detecting BrainVisa Output =============== %

%% ========================== Selecting ROIs  ======================= %
switch cad2Save
    case 'auto'
        switch lower(HemiChar)
            case 'l'
                if exist(opts.bvsulcigroup,'file')
                    [groupsIds, groupsNames, correctNodes, labelsIds] = Detecting_BrainVisa_Groups(opts.bvsulcigroup, 'left');
                    correctNodes = correctNodes(1:length(labelsIds),:);
                else
                    correctNodes = unique(SulcLabels,'rows');
                    correctNodes(find(ismember(correctNodes(:,1:7),'unknown','rows')),:) = [];
                end
                hemiHeader = 'Left-Hemisphere';
            case 'r'
                if exist(opts.bvsulcigroup,'file')
                    [groupsIds, groupsNames, correctNodes, labelsIds] = Detecting_BrainVisa_Groups(opts.bvsulcigroup, 'right');
                    correctNodes = correctNodes(1:length(labelsIds),:);
                else
                    correctNodes = unique(SulcLabels,'rows');
                    correctNodes(find(ismember(correctNodes(:,1:7),'unknown','rows')),:) = [];
                end
                hemiHeader = 'Right-Hemisphere';
        end
        colHeaders = correctNodes;
    case 'lobar'
        switch lower(HemiChar)
            case 'l'
                if exist(opts.bvsulcigroup,'file')
                    [correctNodes, colHeaders] = Detecting_BrainVisa_Groups(opts.bvsulcigroup, 'left');
                else
                    correctNodes = unique(SulcLabels,'rows');
                    correctNodes(find(ismember(correctNodes(:,1:7),'unknown','rows')),:) = [];
                end
                hemiHeader = 'Left-Hemisphere';
            case 'r'
                if exist(opts.bvsulcigroup,'file')
                    [correctNodes, colHeaders] = Detecting_BrainVisa_Groups(opts.bvsulcigroup, 'right');
                else
                    correctNodes = unique(SulcLabels,'rows');
                    correctNodes(find(ismember(correctNodes(:,1:7),'unknown','rows')),:) = [];
                end
                hemiHeader = 'Right-Hemisphere';
        end
    otherwise
        correctNodes = unique(SulcLabels,'rows');
        ind2rem = ismember(correctNodes(:,1:8),'untitled','rows');
        correctNodes(ind2rem,:) = [];
        colHeaders = correctNodes;
        switch lower(HemiChar)
            case 'l'
                hemiHeader = 'Left-Hemisphere';
            case 'r'
                hemiHeader = 'Right-Hemisphere';
        end
end

%% ====================== End of Selecting ROIs  ==================== %

%% ====================== Processing each sulcus =================== %%
[Nstruct, Ncols] = size(correctNodes);
[NrowsLab, NcolsLab] = size(SulcLabels);
if NcolsLab < Ncols
    SulcLabels = [SulcLabels repmat(' ',[NrowsLab Ncols-NcolsLab])];
end

Nsulc   = size(correctNodes,1);
cont = 0;
for st = 1:Nsulc
    Slabelo = deblank(correctNodes(st,:));
    if st == 60
        a = 1;
    end
    % --- Creating varnames
    stVarName = deblank(colHeaders(st,:)); % Structure name
    if st == 1
        measNames = strvcat(['geodDepthMax-' stVarName  '(mm)'],...
            ['geodDepthMin-' stVarName  '(mm)'],...
            ['geodDepthMean-'  stVarName  '(mm)'],...
            ['topLength-'  stVarName  '(mm)'],...
            ['bottomLength-' stVarName  '(mm)'],...
            ['intercepLength-' stVarName  '(mm)'],...
            ['lengthMean-' stVarName  '(mm)'],...
            ['spamMax-' stVarName  '(mm)'],...
            ['spamMin-' stVarName  '(mm)'],...
            ['spamMean-' stVarName  '(mm)']);
    else
        measNames = strvcat(measNames, ['geodDepthMax-' stVarName  '(mm)'],...
            ['geodDepthMin-' stVarName  '(mm)'],...
            ['geodDepthMean-'  stVarName  '(mm)'],...
            ['topLength-'  stVarName  '(mm)'],...
            ['bottomLength-' stVarName  '(mm)'],...
            ['intercepLength-' stVarName  '(mm)'],...
            ['lengthMean-' stVarName  '(mm)'],...
            ['spamMax-' stVarName  '(mm)'],...
            ['spamMin-' stVarName  '(mm)'],...
            ['spamMean-' stVarName  '(mm)']);
    end
    
    % Selecting Sulcus Nodes
    L = length(Slabelo);
    inds = find(ismember(SulcLabels(:,1:L),Slabelo,'rows'));
    NsulcMetrics.length.toplength = 0;
    NsulcMetrics.length.botlength = 0;
    NsulcMetrics.length.meanlength = 0;
    NsulcMetrics.length.intlength = 0;
    NsulcMetrics.depth.profile = 0;
    NsulcMetrics.width.profile = 0;
    if ~isempty(inds)
        
        Nn = length(inds); % Number of Nodes in the sulci mesh
        for indsulc = 1:Nn % Processing each node
            Sulcmet = sulcMetrics{TmtktriIds(inds(indsulc))};
            if ~isempty(sulciLines{TmtktriIds(inds(indsulc))})
                NsulcMetrics.length.toplength = NsulcMetrics.length.toplength + Sulcmet.length.profile(1,end);
                NsulcMetrics.length.botlength = NsulcMetrics.length.botlength + Sulcmet.length.profile(end,end);
                NsulcMetrics.length.meanlength = NsulcMetrics.length.meanlength + mean(Sulcmet.length.profile(:,end));
                NsulcMetrics.length.intlength = NsulcMetrics.length.intlength + Sulcmet.length.measures(1);
                
                % Joining Depth Measures
                NsulcMetrics.depth.profile = [NsulcMetrics.depth.profile;Sulcmet.depth.profile(end,:)'];
                
                % Joining Width Measures
                NsulcMetrics.width.profile = [NsulcMetrics.width.profile;max(Sulcmet.width.profile)'];
            end
        end
    end
    indrem = find(NsulcMetrics.width.profile == 0);
    NsulcMetrics.width.profile(indrem) = [];
    
    indrem = find(NsulcMetrics.depth.profile == 0);
    NsulcMetrics.depth.profile(indrem) = [];
    if  ~isempty(NsulcMetrics.width.profile)
        NsulcMetrics.depth.measures = [min(NsulcMetrics.depth.profile); max(NsulcMetrics.depth.profile);  mode(NsulcMetrics.depth.profile) ; median(NsulcMetrics.depth.profile) ;  mean(NsulcMetrics.depth.profile);  std(NsulcMetrics.depth.profile)];
        NsulcMetrics.width.measures = [min(NsulcMetrics.width.profile); max(NsulcMetrics.width.profile);  mode(NsulcMetrics.width.profile) ; median(NsulcMetrics.width.profile) ;  mean(NsulcMetrics.width.profile);  std(NsulcMetrics.width.profile)];
    else
        NsulcMetrics.width.measures = zeros(1,6);
        NsulcMetrics.width.profile = 0;
        
        % Sulcal Length
        NsulcMetrics.length.measures = zeros(1,7);
        NsulcMetrics.length.profile = 0;
        
        % Sulcal Depth
        NsulcMetrics.depth.measures = zeros(1,6);
        NsulcMetrics.depth.profile = 0;
    end
    
    TempName = [OutdirSPAM filesep subjId '_' HemiChar '_' Slabelo '_' cad2Save '_width.txt'];
    fid  = fopen(TempName,'wt');
    
    line2write = [strvcat('Sulcal_Label',Slabelo) repmat('     ',[2 1 ]) ...
        strvcat('Number_of_Points',num2str(length(NsulcMetrics.width.profile))) repmat('     ',[2 1 ]) ...
        strvcat('Min_Sulcal_Width(mm)',num2str(NsulcMetrics.width.measures(1)))  repmat('     ',[2 1 ]) ...
        strvcat('Max_Sulcal_Width(mm)',num2str(NsulcMetrics.width.measures(2)))  repmat('     ',[2 1 ]) ...
        strvcat('Mode_Sulcal_Width(mm)',num2str(NsulcMetrics.width.measures(3)))  repmat('     ',[2 1 ]) ...
        strvcat('Median_Sulcal_Width(mm)',num2str(NsulcMetrics.width.measures(4)))  repmat('     ',[2 1 ]) ...
        strvcat('Mean_Sulcal_Width(mm)',num2str(NsulcMetrics.width.measures(5)))  repmat('     ',[2 1 ]) ...
        strvcat('Std_Sulcal_Width(mm)',num2str(NsulcMetrics.width.measures(6)))];
    
    fprintf(fid,'%s\n', line2write(1,:));
    fprintf(fid,'%s\n', line2write(2,:));
    fclose(fid);
    
    % Saving Depth Results
    TempName = [OutdirDEPTH filesep subjId '_' HemiChar '_' Slabelo '_' cad2Save '_depth.txt'];
    fid  = fopen(TempName,'wt');
    
    line2write = [strvcat('Sulcal_Label',Slabelo) repmat('     ',[2 1 ]) ...
        strvcat('Number_of_DepthProfiles',num2str(length(NsulcMetrics.depth.profile))) repmat('     ',[2 1 ]) ...
        strvcat('Min_Sulcal_Depth(mm)',num2str(NsulcMetrics.depth.measures(1)))  repmat('     ',[2 1 ]) ...
        strvcat('Max_Sulcal_Depth(mm)',num2str(NsulcMetrics.depth.measures(2)))  repmat('     ',[2 1 ]) ...
        strvcat('Mode_Sulcal_Depth(mm)',num2str(NsulcMetrics.depth.measures(3)))  repmat('     ',[2 1 ]) ...
        strvcat('Median_Sulcal_Depth(mm)',num2str(NsulcMetrics.depth.measures(4)))  repmat('     ',[2 1 ]) ...
        strvcat('Mean_Sulcal_Depth(mm)',num2str(NsulcMetrics.depth.measures(5)))  repmat('     ',[2 1 ]) ...
        strvcat('Std_Sulcal_Depth(mm)',num2str(NsulcMetrics.depth.measures(6)))];
    fprintf(fid,'%s\n', line2write(1,:));
    fprintf(fid,'%s\n', line2write(2,:));
    fclose(fid);
    
    % Saving Length Results
    TempName = [OutdirLENGTH filesep subjId '_' HemiChar '_' Slabelo '_' cad2Save '_length.txt'];
    fid  = fopen(TempName,'wt');
    
    line2write = [strvcat('Sulcal_Label',Slabelo) repmat('     ',[2 1 ]) ...
        strvcat('Top_Length(mm)',num2str(NsulcMetrics.length.toplength)) repmat('     ',[2 1 ]) ...
        strvcat('Bottom_Length(mm)',num2str(NsulcMetrics.length.botlength))  repmat('     ',[2 1 ]) ...
        strvcat('Intercept_Length(mm)',num2str(NsulcMetrics.length.intlength))  repmat('     ',[2 1 ]) ...
        strvcat('Mean_Length(mm)',num2str(NsulcMetrics.length.meanlength))  repmat('     ',[2 1 ])];
    
    fprintf(fid,'%s\n', line2write(1,:));
    fprintf(fid,'%s\n', line2write(2,:));
    fclose(fid);
    %     else
    %         cont = cont + 1;
    %     end
    if st == 1
        metricsMat = [NsulcMetrics.depth.measures(2) NsulcMetrics.depth.measures(1) NsulcMetrics.depth.measures(5) ...
            NsulcMetrics.length.toplength NsulcMetrics.length.botlength NsulcMetrics.length.intlength NsulcMetrics.length.meanlength ...
            NsulcMetrics.width.measures(2) NsulcMetrics.width.measures(1) NsulcMetrics.width.measures(5)];
    else
        metricsMat = [metricsMat NsulcMetrics.depth.measures(2) NsulcMetrics.depth.measures(1) NsulcMetrics.depth.measures(5) ...
            NsulcMetrics.length.toplength NsulcMetrics.length.botlength NsulcMetrics.length.intlength NsulcMetrics.length.meanlength ...
            NsulcMetrics.width.measures(2) NsulcMetrics.width.measures(1) NsulcMetrics.width.measures(5)];
    end
end

switch cad2Save
    % Grouping into lobes
    case 'auto'
        % --------------- Grouping into Lobes
        sts = unique(labelsIds);
        Nlobes = length(sts);
        for j = 1:Nlobes
            ind = find(labelsIds == sts(j));
            
            % Depth
            T                = nonzeros(metricsMat(:,(ind-1)*10+1));
            if T==0
                LmaxDepth        = 0;
                LminDepth        = 0;
                LmeanDepth       = 0;
                
                LtopLength       = 0;
                LbottomLength    = 0;
                LinterLength     = 0;
                LmeanLength      = 0;

                LmaxSPAM         = 0;
                LminSPAM         = 0;
                LmeanSPAM        = 0;
            else
                LmaxDepth        = max(T);
                T                = nonzeros(metricsMat(:,(ind-1)*10+2));
                LminDepth        = min(T);
                T                = metricsMat(:,(ind-1)*10+3);
                LmeanDepth       = sum(T,2)./(sum(logical(T),2) + eps);
                
                % Length
                LtopLength       = sum(metricsMat(:,(ind-1)*10+4),2);
                LbottomLength    = sum(metricsMat(:,(ind-1)*10+5),2);
                LinterLength     = sum(metricsMat(:,(ind-1)*10+6),2);
                LmeanLength      = sum(metricsMat(:,(ind-1)*10+7),2);
                
                % Width
                T                = nonzeros(metricsMat(:,(ind-1)*10+8));
                LmaxSPAM         = max(T);
                T                = nonzeros(metricsMat(:,(ind-1)*10+9));
                LminSPAM         = min(T);
                T                = metricsMat(:,(ind-1)*10+10);
                LmeanSPAM        = sum(T,2)./(sum(logical(T),2) + eps);
            end
            if j == 1
                tempMatL = [LmaxDepth LminDepth LmeanDepth LtopLength LbottomLength LinterLength LmeanLength LmaxSPAM LminSPAM LmeanSPAM];
                measNamesL = strvcat(['geodDepthMax-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['geodDepthMin-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['geodDepthMean-'  deblank(groupsNames(j,:))  '(mm)'],...
                    ['topLength-'  deblank(groupsNames(j,:))  '(mm)'],...
                    ['bottomLength-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['intercepLength-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['lengthMean-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['spamMax-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['spamMin-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['spamMean-' deblank(groupsNames(j,:))  '(mm)']);
            else
                tempMatL = [tempMatL LmaxDepth LminDepth LmeanDepth LtopLength LbottomLength LinterLength LmeanLength LmaxSPAM LminSPAM LmeanSPAM];
                measNamesL = strvcat(measNamesL, ['geodDepthMax-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['geodDepthMin-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['geodDepthMean-'  deblank(groupsNames(j,:))  '(mm)'],...
                    ['topLength-'  deblank(groupsNames(j,:))  '(mm)'],...
                    ['bottomLength-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['intercepLength-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['lengthMean-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['spamMax-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['spamMin-' deblank(groupsNames(j,:))  '(mm)'],...
                    ['spamMean-' deblank(groupsNames(j,:))  '(mm)']);
            end
        end
        measNames = strvcat(measNamesL,measNames);
        metricsMat = [tempMatL metricsMat];
    case 'lobar'
        tempMatL = metricsMat;
end

% Creating Hemisphere Values
measNames = strvcat(['geodDepthMax-' hemiHeader  '(mm)'],...
    ['geodDepthMin-' hemiHeader  '(mm)'],...
    ['geodDepthMean-'  hemiHeader  '(mm)'],...
    ['topLength-'  hemiHeader  '(mm)'],...
    ['bottomLength-' hemiHeader  '(mm)'],...
    ['intercepLength-' hemiHeader  '(mm)'],...
    ['lengthMean-' hemiHeader  '(mm)'],...
    ['spamMax-' hemiHeader  '(mm)'],...
    ['spamMin-' hemiHeader  '(mm)'],...
    ['spamMean-' hemiHeader  '(mm)'],...
    measNames);
metricsMat = [max(tempMatL(1:10:end)) ...
    min(tempMatL(2:10:end)) ...
    mean(tempMatL(3:10:end)) ...
    sum(tempMatL(4:10:end)) ...
    sum(tempMatL(5:10:end)) ...
    sum(tempMatL(6:10:end)) ...
    sum(tempMatL(7:10:end)) ...
    max(tempMatL(8:10:end)) ...
    min(tempMatL(9:10:end))...
    mean(tempMatL(10:10:end)) ...
    metricsMat];

%% ==================== End Processing each sulcus  ==================== %%
varargout{1} = metricsMat;
varargout{2} = measNames;
return;
