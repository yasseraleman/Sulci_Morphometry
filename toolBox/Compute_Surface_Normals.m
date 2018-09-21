function varargout = Compute_Surface_Normals(varargin);
%
% Syntax :
%     Surfout = Compute_Surface_Normals(Surf);
%
% This function computes surface normals.
%
% Input Parameters:
%        Surf                   : Surface variable (file, struct or cellarray).
%
% Output Parameters:
%        Surfout                : Output Surface variable.
%
% See also: Plot_Surf Surf_Color Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%
if nargin == 0
   error('One Input is mandatory');
   return; 
end
if nargin > 1
    error('To Many Input Parameters');
    return;
end
if nargout > 1
    error('To Many Output Parameters');
    return;
end
Surf = varargin{1};
%% ========================= End of Checking Inputs ==================== %%

%% ============================= Main Program ========================== %%
if ischar(Surf) % Verifying if the Surf variable is a text
    if exist(Surf,'file') % Verifying if the surface file exist
        try
            [OutFiles, SurfF] = Exp_Surf(Surf, '0', '','', 'imp','n'); % Reading surface
            Surfout = SurfF{1};
        catch
            error('Unrecognized Surface format');
            return;
        end
    end
elseif isstruct(Surf); % Verifying if the surface is a matlab structure
    if isfield(Surf,'SurfData')% Verifying important fields
        if ~isfield(Surf.SurfData,'faces')|~isfield(Surf.SurfData,'vertices'); % Verifying important fields
            error('Important fields from the Surface matlab are Missing');
            return
        else
            Surfout = Surf;
        end
    else
        error('Unrecognized Surface format');
        return;
    end
elseif iscell(Surf); % Verifying if the surface is a matlab cell
    for i = 1:length(Surf)
        Surfttemp = Surf{i};
        if isfield(Surfttemp,'SurfData')% Verifying important fields
            if ~isfield(Surfttemp.SurfData,'faces')|~isfield(Surfttemp.SurfData,'vertices'); % Verifying important fields
                error('Important fields from the Surface matlab are Missing');
                return
            else
                Surfout{i,1} = Surfttemp;
            end
        else
            error('Unrecognized Surface format');
            return;
        end
    end
else
    error('Unrecognized Surface format');
    return;
end


% Computing Normals
Surfout.SurfData.VertexNormals = patchnormals(Surfout.SurfData);

in = Check_Points_Inside_Surface(Surfout,Surfout.SurfData.vertices(1,:) + Surfout.SurfData.VertexNormals(1,:));
if in == 1; % Vertices outside hull surface
    Surfout.SurfData.VertexNormals = Surfout.SurfData.VertexNormals*-1; % Flip normals
end


%% ========================= End of Main Program ======================= %%
% Outputs
varargout{1} = Surfout;
return