function varargout = Surface_Checking(varargin);
%
% Syntax :
%     Surfout = Surface_Checking(Surf);
%
% This function checks surface format before reading it.
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
Surfout = '';
Surf = varargin{1};
%% ========================= End of Checking Inputs ==================== %%

%% ============================= Main Program ========================== %%
if ischar(Surf) % Verifying if the Surf variable is a text
    if exist(Surf,'file') % Verifying if the surface file exist
        try
            [OutFiles, SurfF] = Exp_Surf(Surf, '0', '','', 'imp','n'); % Reading surface
            Surfout = SurfF{1};
            if ~isfield(Surfout.SurfData,'VertexNormals')&exist(which('Compute_Surface_Normals.m'),'file')
                Surfout = Compute_Surface_Normals(Surfout);
            end
        catch
            error('Unrecognized Surface format');
            return;
        end
    end
elseif isstruct(Surf); % Verifying if the surface is a matlab structure
   if isfield(Surf(1),'SurfData')% Verifying important fields
        if ~isfield(Surf(1).SurfData,'faces')|~isfield(Surf(1).SurfData,'vertices'); % Verifying important fields
            error('Important fields from the Surface matlab are Missing');
            return
        else
            if ~isfield(Surf(1).SurfData,'VertexNormals')&exist(which('Compute_Surface_Normals.m'),'file')
                if size(Surf(1).SurfData.faces,2) ==3
                    Surf(1) = Compute_Surface_Normals(Surf(1));
                end
            end
            Surfout = Surf;
        end
    else
        error('Unrecognized Surface format');
        return;
    end
elseif iscell(Surf); % Verifying if the surface is a matlab cell
     Surf = Surf(:);
     indempty = cellfun(@isempty,Surf); % Remove empty cells
     Surf(indempty) = [];
    for i = 1:length(Surf)
        Surfttemp = Surf{i};
        if isfield(Surfttemp(1),'SurfData')% Verifying important fields
            if ~isfield(Surfttemp(1).SurfData,'faces')|~isfield(Surfttemp(1).SurfData,'vertices'); % Verifying important fields
                error('Important fields from the Surface matlab are Missing');
                return
            else
                if ~isfield(Surfttemp(1).SurfData,'VertexNormals')&exist(which('Compute_Surface_Normals.m'),'file')
                    if size(Surfttemp(1).SurfData.faces,2) ==3
                        Surfttemp(1) = Compute_Surface_Normals(Surfttemp(1));
                    end
                end
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
%% ========================= End of Main Program ======================= %%
% Outputs
varargout{1} = Surfout;
return