function varargout = Read_Surface(varargin);
%
% Syntax :
%  Surf = Read_Surface(surfFile, charVal);
%
% This script reads any surface format and a specified overlay.
%
% Input Parameters:
%       surfFile              : Surface Filename
%
%       charVal               : Map Filename
%
% Output Parameters:
%      Surf                   : Surface struct
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% July 22th 2014
% Version $1.0

%% ===================== Checking Input Parameters ====================== %
if nargin > 1
    surfFile = varargin{1};
    charVal = varargin{2};
elseif nargin == 1
    surfFile = varargin{1};
else
    error('Too many inputs');
    return
end
%% ===================== End of Checking Input Parameters =============== %

%% ========================== Main Program ============================== %

% ---------------------- Reading Surface File --------------------------- %
try
    [OutFiles, SurfF] = Exp_Surf(surfFile, '0', '','', 'imp','n');
    Surf = SurfF{1};
catch
    
    error('Unrecognized Surface File');
end

% ------------------ Reading Characteristic Map ------------------------- %
if nargin > 1
    try
        
        if ischar(charVal)
            if exist(charVal,'file')
                [mapValues, ctab] = read_cfiles(charVal);
                boolVar = 1;
            end
            if isempty(strfind(charVal,filesep))
                ind = strfind(surfFile,'.');
                if ~isempty(ind)
                    switch charVal
                        case 'thickness'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'curv'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'pial_lgi'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'aparc.annot'
                            [pth, nm, ext] = fileparts(surfFile);
                            [tempName, nmt, ext] = fileparts(pth);
                            if length(nm) >1
                                hemi = nm(1:2);
                            else
                                hemi = nm;
                            end
                            mapTemp = [tempName filesep 'label' filesep hemi '.' charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'aparc.a2009s.annot'
                            [pth, nm, ext] = fileparts(surfFile);
                            [tempName, nmt, ext] = fileparts(pth);
                            if length(nm) >1
                                hemi = nm(1:2);
                            else
                                hemi = nm;
                            end
                            mapTemp = [tempName filesep 'label' filesep hemi '.' charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'sulc'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'area'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'volume'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'myelin'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'fa'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'md'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'vf'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'cp'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'cs'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'ga'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'ad'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'rd'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        case 'ra'
                            mapTemp = [surfFile(1:ind) charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                        otherwise
                            [pth, nm, ext] = fileparts(surfFile);
                            [tempName, nmt, ext] = fileparts(pth);
                            if length(nm) >1
                                hemi = nm(1:2);
                            else
                                hemi = nm;
                            end
                            mapTemp = [tempName filesep 'label' filesep hemi '.' charVal];
                            if exist(mapTemp,'file')
                                [mapValues, ctab] = read_cfiles(mapTemp);
                                boolVar = 1;
                            else
                                boolVar = 0;
                            end
                    end
                end
            end
        elseif isnumeric(charVal)
            mapValues = charVal;
            ctab.table = 0;
        end
        if size(Surf.SurfData.vertices,1) == size(mapValues(:),1)
            if boolVar
                Surf.Is = mapValues;
                if sum(ctab.table(:)) ~=0 % Just for FreeSurfer .annot files
                    if isfield(ctab,'struct_names')
                        
                        % Removing Medial Wall vertices
%                         tempname = char(ctab.struct_names);
%                         indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
%                         indcc = find(ismember(tempname(:,1:14),'corpuscallosum','rows') == 1);
%                         ctab.table([indu indcc],:) = [];
%                         ctab.struct_names([indu indcc]) = [];
%                         mwallind = find(ismember(mapValues,ctab.table(:,5)) == 0); % Medial Wall Indexes
                        
                        Colors = ctab.table(:,1:3); % FreeSurfer Colors extrated from .annot file
                        
                        sts = ctab.table(:,5);   % Unifying structures labels
                        
                        Surf.SurfData.FaceVertexCData = .41*ones(length(mapValues),3); % Gray Values to undefined regions
                        for i = 1:length(sts)
                            inds = find(mapValues == sts(i));
                            Surf.SurfData.FaceVertexCData(inds,:) = repmat(Colors(i,:)/255, [length(inds) 1]);
                        end
                        %Surf.SurfData.FaceVertexCData(mwallind,:) = ones(length(mwallind),3);
                    end
                end
            end
        else
            disp('WARNING: Map length is different from the number of vertices');
        end
    catch
        error('Unrecognized Map File');
    end
end
varargout{1} = Surf; % Output Surface
%% ====================== End of Main Program =========================== %
return;