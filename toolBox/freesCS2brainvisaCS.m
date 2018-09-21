function varargout = freesCS2brainvisaCS(varargin);
%
% Syntax :
% Surf = freesCS2brainvisaCS(Surf,Imfile,ttype);
%
% This script move Surfaces in freesurfer coordinate system to Brainvisa 
% coordinate system.
%
% Input Parameters:
%   Surf   : Individual surface structure. Surface coordinates must be in
%             voxels
%   Imfile  : Image file in FreeSurfer space
%   ttype   : Transformation type(ie f2b: freesurfer to brainvisa, 
%                                    b2f: brainvisa to freesurfer)
%  
% Output Parameters:
%   Surfi    : Transformed Individual surface structure
%
% Related references:
%
% See also:  save_texBrainvisa
%__________________________________________________
% Authors:  Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% April 16th 2012
% Version $1.0

warning off;
if nargin <2
    error(' There are missing inputs');
    return;
end
Surf = varargin{1};
Imfile = varargin{2};
if nargin <3
    ttype = 'f2b';
else
    ttype = varargin{3};
end

if ischar(Imfile)
    V = spm_vol(Imfile);
elseif isstruct(Imfile)
    if isfield(Imfile,'mat')
        V = Imfile;
    else
        error('Invalid Image File or volume struct');
        return;
    end
else
    error('Invalid Image File or volume struct');
    return;
end

Mat = [1 0 0 -1;0 0 -1 V(1).dim(2);0 1 0 -1;0 0 0 1];
Ns = length(Surf);
switch ttype
    case 'f2b'
        for i = 1:Ns
            vertvox = (Mat*inv(V.mat)*[Surf(i).SurfData.vertices ones(size(Surf(i).SurfData.vertices,1),1)]')';
            Surf(i).SurfData.vertices = vertvox(:,1:3);
        end
    case 'b2f'
        for i = 1:Ns
            vertvox = (V.mat*inv(Mat)*[Surf(i).SurfData.vertices ones(size(Surf(i).SurfData.vertices,1),1)]')';
            Surf(i).SurfData.vertices = vertvox(:,1:3);
        end
end
varargout{1} = Surf;
return;