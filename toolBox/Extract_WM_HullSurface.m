function [varargout] = Extract_WM_HullSurface(varargin);
%
% Syntax :
%   [HullSurfaces] = Extract_WM_HullSurface(FreeSurferDatabaseDir,Idfile);
%
% This script creates hemispheric white matter hull surface from freesurfer parcellation
% atlas. 
%
% Input Parameters:
%   FreeSurferDatabaseDir   :  FreeSurfer Output directory.
%     Idfile                :  Ids File
%
% Output Parameters:
%     HullSurfaces          :  White Matter hull surfaces filename.
%
% Related references:
%
%
% See also:
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% May 12th 2012
% Version $1.0


% Idfile = '1Test_HCP_899885-20140807-T1wMPR1';
% FreeSurferDatabaseDir = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing';

%% =================== Checking Input Parameters ======================== %

if nargin == 2
    FreeSurferDatabaseDir = varargin{1};
    Idfile = varargin{2};
else
    error('Wrong Inputs');
    return;
end
%% =================== End of Checking Input Parameters ================= %
Ndil = 5;
if exist(Idfile,'file')
    Ids = char(textread(Idfile,'%s'));
else
    Ids = Idfile;
end
Ns = size(Ids,1);
HullSurfaces = '';
for i = 1:Ns
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    
    % ---- Creating Output Directories
    
    %% ====================== General Processing ======================= %%
    lhullfilename = [FreeSurferDatabaseDir filesep subjId filesep 'surf' filesep 'lh.wmhull'];
    ImageR = [FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'ribbon.mgz'];
    
    % Converting FreeSurfer Ribbon Image to nifti
    system(['rm -rf ' FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'ribbon.nii']);
    system(['rm -rf ' FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'fribbon.nii']);
    cad = ['mri_convert -i ' FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'ribbon.mgz -o ' FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'ribbon.nii'] ;
    system(cad);
    Vfs = spm_vol([FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'ribbon.nii']);
    outputs = flip_to_axial_nii(Vfs.fname,0); % Flipping to axial orientation
    outputs = [FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'fribbon.nii'];
    V = spm_vol(outputs);
    I = int16(spm_read_vols(V));
    
    %% ========== Creating Left hemisphere wm hull surface ============= %%
    ind = find(I == 2); % Left hemisphere white matter voxels
    It = logical(zeros(size(I)));
    It(ind) = 1;  %
    
    Id = logical(zeros(size(I))); % Image to dilate
    ind = find((I == 42)|(I == 41)); %  Right hemisphere voxels (GM+ WM)
    Id(ind) = 1;
    for i = 1:Ndil
        Id = imdilate(Id,strel(ones(3,3,3)));
    end
    
    ind2rem= find(Id.*It);
    It(ind2rem) = 0; % Remove voxels from the interface between hemispheres
    
    %% ============== White Matter Hull Extraction ===================== %%
    %  Extracting Hull
    
    It = Iso_Rem(It,13);
    HullSurfMat = Extract_Hull_Surface_from_Mask(It);
    
    % Hull vertices in mm
    voxs = (V.mat*[HullSurfMat.SurfData.vertices(:,1) HullSurfMat.SurfData.vertices(:,2) HullSurfMat.SurfData.vertices(:,3) ones(size(HullSurfMat.SurfData.vertices(:,1),1),1)]')';
    HullSurfMat.SurfData.vertices = [voxs(:,1) voxs(:,2) voxs(:,3)];
    HullSurfMat.SurfData.faces = HullSurfMat.SurfData.faces;
    
    % Hull Smoothing
% %     FV=smoothpatch(HullSurfMat.SurfData,1,50);
% %     HullSurfMat.SurfData.vertices = FV.vertices;
    

    
    HullSurfMat = Compute_Surface_Normals(HullSurfMat);
    
    distvert = sqrt(sum(sqrt(sum(V.mat(1:3,1:3).^2)).^2))/4;
    HullSurfMat.SurfData.vertices = HullSurfMat.SurfData.vertices - HullSurfMat.SurfData.VertexNormals*distvert;
    OutFile = save_free_surf(HullSurfMat, lhullfilename);
    
    
    %% ============== End of  White Matter Hull Extraction ============= %%
    HullSurfaces = strvcat(HullSurfaces,lhullfilename);
    
    %% ======== End of Creating Left hemisphere wm hull surface ======== %%
    
    
    %% ============ Creating Right hemisphere wm hull surface ========== %%
    rhullfilename = [FreeSurferDatabaseDir filesep subjId filesep 'surf' filesep 'rh.wmhull'];
    ind = find(I == 41); % Left hemisphere white matter voxels
    It = logical(zeros(size(I)));
    It(ind) = 1;  %
    
    Id = logical(zeros(size(I))); % Image to dilate
    ind = find((I == 2)|(I == 3)); %  Right hemisphere voxels (GM+ WM)
    Id(ind) = 1;
    for i = 1:Ndil
        Id = imdilate(Id,strel(ones(3,3,3)));
    end
    
    ind2rem= find(Id.*It);
    It(ind2rem) = 0; % Remove voxels from the interface between hemispheres
    
    %% ============== White Matter Hull Extraction ===================== %%
    %  Extracting Hull
    HullSurfMat = Extract_Hull_Surface_from_Mask(It);
    
    % Hull vertices in mm
    voxs = (V.mat*[HullSurfMat.SurfData.vertices(:,1) HullSurfMat.SurfData.vertices(:,2) HullSurfMat.SurfData.vertices(:,3) ones(size(HullSurfMat.SurfData.vertices(:,1),1),1)]')';
    HullSurfMat.SurfData.vertices = [voxs(:,1) voxs(:,2) voxs(:,3)];
    HullSurfMat.SurfData.faces = HullSurfMat.SurfData.faces;
    
    % Hull Smoothing
% %     FV=smoothpatch(HullSurfMat.SurfData,1,50);
% %     HullSurfMat.SurfData.vertices = FV.vertices;
    
   OutFile = save_free_surf(HullSurfMat, rhullfilename);
    %% ============== End of White Matter Hull Extraction ============= %%
    
    %% ======== End of Creating Right hemisphere wm hull surface ======== %%
    HullSurfaces = strvcat(HullSurfaces,rhullfilename);
end
varargout{1} = HullSurfaces;
return;