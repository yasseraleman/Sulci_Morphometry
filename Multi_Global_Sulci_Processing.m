function [varargout] = Multi_Global_Sulci_Processing(varargin);
%
% Syntax :
%  [statFiles] = Multi_Global_Sulci_Processing(Ids, opts);
%
% This script computes sulci metrics for different subjects
%
% Input Parameters:
%       IdFile                  : Ids Filename. It also can be a list of
%                                 subjects Ids specified by a string list of Ids.
%
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
% See also: Whole_Brain_Sulci_Processing Global_Sulci_Processing Hemispheric_Sulci_Processing Sulci_Nodal_Processing  Compute_Node_Metrics
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0


% % % % % %% ==================== Checking Input Parameters ======================= %
% % % % % if nargin <1
% % % % %     error('One input is mandatory');
% % % % %     return;
% % % % % end
% % % % % IdFile = varargin{1};
% % % % % if nargin == 2
% % % % %     opts = varargin{2};
% % % % % end
% % % % % if ~isfield(opts,'brainvisadir')     %(default BrainVisa _*Hemi.gii)
% % % % %     error('BrainVisa Directory is mandatory');
% % % % %     return;
% % % % % else
% % % % %     if ~exist(opts.brainvisadir, 'dir')
% % % % %         error('BrainVisa Directory does not exist');
% % % % %         return;
% % % % %     end
% % % % % end
% % % % % if ~isfield(opts,'subjid')
% % % % %     error('Subject ID is mandatory');
% % % % %     return;
% % % % % else
% % % % %     if ~exist([opts.brainvisadir filesep 'subjects' filesep opts.subjid], 'dir')
% % % % %         error('Subject Id does not exist inside the BrainVisa Directory');
% % % % %         return;
% % % % %     end
% % % % % end
% % % % % 
% % % % % if ~isfield(opts,'pialsurf')
% % % % %     opts.pialsurf = 'brainvisa'; % Pial Surface file (default BrainVisa _*Hemi.gii)
% % % % % end
% % % % % if ~isfield(opts,'hullsurf')     %(default BrainVisa _*Hemi.gii)
% % % % %     opts.hullsurf = 'brainvisa';
% % % % % end
% % % % % if ~isfield(opts,'pialsurf')
% % % % %     opts.pialsurf = 'brainvisa'; % Pial Surface file (default BrainVisa _*Hemi.gii)
% % % % % end
% % % % % if ~isfield(opts,'hemisphere') % Hemisphere
% % % % %     opts.hemisphere = 'left';
% % % % %     warning('Left Hemisphere has been selected as default');
% % % % % end
% % % % % if ~isfield(opts,'maxsulcawidth') % mm. Maximum Sulcal Width
% % % % %     opts.maxsulcawidth = 12;
% % % % % end
% % % % % if ~isfield(opts,'mindist2hull') % mm. Minimum distance to the hull surface mm
% % % % %     opts.mindist2hull = 5;
% % % % % end
% % % % % if ~isfield(opts,'mindist2hull') % degrees. Minimum angle allowed between normals
% % % % %     opts.mindist2hull = 5;
% % % % % end
% % % % % if ~isfield(opts,'ncurvp') % Number of curve points
% % % % %     opts.ncurvp = 40;
% % % % % end
% % % % % 
% % % % % if strcmp(opts.hullsurf,'freesurfer')
% % % % %     if ~isfield(opts,'freesurferdir')     %(default BrainVisa _*Hemi.gii)
% % % % %         error('Freesurferdir Directory is mandatory');
% % % % %         return;
% % % % %     else
% % % % %         if ~exist(opts.freesurferdir, 'dir')
% % % % %             error('Freesurferdir Directory does not exist');
% % % % %             return;
% % % % %         end
% % % % %     end
% % % % %     if ~exist([opts.freesurferdir filesep opts.subjid], 'dir')
% % % % %         error('Subject Id does not exist inside the FreeSurfer Directory');
% % % % %         return;
% % % % %     end
% % % % % end
% % % % % 
% % % % % %% ==================== End of Checking Input Parameters ================ %
% IdFile = varargin{1};
opts.brainvisadir  =  '/media/COSAS/8-BrainVISADataBase-HCP';
opts.freesurferdir = '/media/HCPData/5-freesurfer_processing';
opts.pialsurf = 'freesurfer';
opts.hullsurf = 'freesurfer';

IdFile = '/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';
% IdFile = strvcat('1SUBJECT_100206_TestSubject');

 %IdFile = '/media/HCPData/30-Ids_textFiles/Ids_BrainVisaProcessing.txt';
% IdFile = '/media/Data/PROCESSING_RESULTS/HCP/8-BrainVisaDataBase/Ids.txt';

if exist(deblank(IdFile),'file');
    Ids = char(textread(IdFile,'%s'));
else
    Ids = IdFile;
end
% Ids = strvcat('101006','104820','108121');

Ns = size(Ids, 1);
statFiles = '';
for i = 202:374
    opts.subjid = deblank(Ids(i,:));
%     lockFile = [opts.brainvisadir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep opts.subjid '_sulcMorph_corr.lock'];
%     if ~exist(lockFile,'file')
%         fid = fopen(lockFile,'wt');
        disp(['Proccessing Subject: ' opts.subjid '  ==> ' num2str(i) ' of ' num2str(Ns)]);
        [statFile] = Whole_Brain_Sulci_Processing(opts);
        statFiles = strvcat(statFiles,statFile);
%         fclose(fid);
%     end
end
varargout{1} = statFiles;
%% %%%%%%%%%%%%%%%%%%%%%%%% End of Statistics %%%%%%%%%%%%%%%%%%%%%% %%
return;