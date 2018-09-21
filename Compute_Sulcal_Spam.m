function varargout = Compute_Sulcal_Spam(varargin)
%
% Syntax :
%       [Surfwidth, sulcalWidth] = Compute_Sulcal_Spam(PialSurfMat, Surfo, opts);
%
% This function computes sulcal spam intercepting the median mesh normals
% in the Pial surface
%
%
% Input Parameters:
%       PialSurfMat             : Pial Surface (Matlab Format)
%       Surfo                   : Median Mesh (Matlab Format).
%
% Output Parameters:
%      Surfwidth                : Sulcal Spam lines. They can be plotted
%                                 using Plot_Surf.
%      sulcalWidth              : Sulcal Width Vector.
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% February 13th 2015
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return
else
    PialSurfMat  = varargin{1};
    Surfo = varargin{2};
    PialSurfMat = Surface_Checking(PialSurfMat);
    Surfo = Surface_Checking(Surfo);
end
%% =========================== Input parameters  =========================%
if nargin < 3
    opts.maxsulcalwidth  = 12; % mm. Maximum Sulcal Width
    opts.epsilon = 10^-5;          % Tolerance to estimate barycentric coordinates
    opts.verbose = 1;              % Maximum Sulcal Width
else
    opts = varargin{3};
    if ~isfield(opts,'maxsulcalwidth')
        opts.maxsulcalwidth = 12;      % Maximum Sulcal Width
    end
    if ~isfield(opts,'epsilon')
        opts.epsilon = 1;              % Tolerance to estimate barycentric coordinates
    end
    if ~isfield(opts,'verbose')
        opts.verbose = 1;              % Maximum Sulcal Width
    end
end
if nargin > 3
    error('To many inputs');
    return;
end
%% ==================== End of Input parameters  =========================%

% % % load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab_test1.mat');
%% ================= Computing Sulcal Spam  ===================== %
if opts.verbose
    disp(' ');
    disp('Computing Width lines ... ');
    tic;
end
cont = 0;

Normals = Surfo.SurfData.VertexNormals;
% Pial Surface planes
PialSurfMat = Compute_Surface_Normals(PialSurfMat);
VertP = PialSurfMat.SurfData.vertices;FacesP = PialSurfMat.SurfData.faces;NormalsP = PialSurfMat.SurfData.VertexNormals;
np = cross(VertP(FacesP(:,1),:)-VertP(FacesP(:,2),:),VertP(FacesP(:,3),:)-VertP(FacesP(:,2),:),2); Dp = -1*dot(np,VertP(FacesP(:,1),:),2);

v0 = VertP(FacesP(:,3),:) - VertP(FacesP(:,1),:);
v1 = VertP(FacesP(:,2),:) - VertP(FacesP(:,1),:);
    
dot00 = dot(v0, v0, 2);
dot01 = dot(v0, v1, 2);
dot11 = dot(v1, v1, 2);

    
P1 = Surfo.SurfData.vertices+Normals; P2 = Surfo.SurfData.vertices-Normals; % Creating perpendicular lines to sulcus walls
Npoints = length(P1);
cont = 0;
Verts = Surfo.SurfData.vertices;

% if Npoints >500
%     P = round(linspace(1,Npoints,500));
%     P1 = P1(P,:);
%     P2 = P2(P,:);
%     Verts = Verts(P,:);
% end
Npoints = length(P1);
sulcalWidth = zeros(Npoints,1);
Surfwidth.SurfData.vertices = zeros(2*Npoints,3);
Surfwidth.SurfData.faces = zeros(Npoints,2);
for po =1:Npoints
    Num = dot(np,repmat(P1(po,:),[size(np,1) 1]),2)+Dp; Den = dot(np,repmat(P2(po,:)-P1(po,:),[size(np,1) 1]),2); % Creating Lines
    
    t = -1*Num./Den;clear Num Den; % Lines parameters
    xint = single(repmat(P1(po,1)',[size(FacesP,1) 1])+t.*(repmat((P2(po,1)-P1(po,1))',[size(FacesP,1) 1]))); % Line parametrization
    yint = single(repmat(P1(po,2)',[size(FacesP,1) 1])+t.*(repmat((P2(po,2)-P1(po,2))',[size(FacesP,1) 1])));
    zint = single(repmat(P1(po,3)',[size(FacesP,1) 1])+t.*(repmat((P2(po,3)-P1(po,3))',[size(FacesP,1) 1])));clear t;
    

    v2 = [xint yint zint]- VertP(FacesP(:,1),:);
    
    % dot products
    dot02 = dot(v0, v2, 2);
    dot12 = dot(v1, v2, 2);
    
    % barycentric coordinates
    invDenom = 1 ./ (dot00 .* dot11 - dot01 .* dot01);
    u = (dot11 .* dot02 - dot01 .* dot12) .* invDenom;
    v = (dot00 .* dot12 - dot01 .* dot02) .* invDenom;
    out = u >= -opts.epsilon & v >= -opts.epsilon & (u+v-opts.epsilon) <= 1;
    ind = find(out);
    
    [inte,indord] = unique([xint(ind) yint(ind) zint(ind)],'rows');ind = ind(indord); % Line interception with pial surface
    
    if ~isempty(inte)
        t = NormalsP(FacesP(ind,:),:);
        intNormals = reshape(mean(reshape(reshape(t,[length(ind) prod(size(t))/length(ind)])',[3 length(ind)*3]))',[3 length(ind)])';
        norma = normm(intNormals);
        intNormals = intNormals./[norma norma norma];
        %     orientsign = sign((dot(inte -repmat(Surf2Process.SurfData.vertices(po,:),[size(inte,1) 1]),repmat(Surf2Process.SurfData.VertexNormals(po,:),[size(inte,1) 1]),2))); % Detecting Orientation
        distances = sqrt(sum((inte - repmat(Verts(po,:),[size(inte,1) 1])).^2,2)); % Computing distance to the pial surface
        
        [valMin,posit] = min(distances);
        MinPositCoord = inte(posit,:);
        t1 = inte - repmat(MinPositCoord,[size(inte,1) 1]);
        t2 = repmat(intNormals(posit,:),[size(inte,1) 1]);
        indPositDir = find(sum(t1.*t2,2) > 0);
        if ~isempty(indPositDir)
            distances = sqrt(sum((inte(indPositDir,:) - repmat(MinPositCoord,[length(indPositDir) 1])).^2,2));
            ind = find((distances == min(nonzeros(distances))));
            MinNegatCoord = inte(indPositDir(ind),:);
            sulcw = distances(ind);
            
            sulcalWidth(po,1) = sulcw;
            %plot3([MinPositCoord(1,1);MinNegatCoord(1,1)],[MinPositCoord(1,2);MinNegatCoord(1,2)],[MinPositCoord(1,3);MinNegatCoord(1,3)],'-','Color',[0 0 1],'Markersize',30);
            Surfwidth.SurfData.vertices(2*po-1:2*po,:) = [MinPositCoord;MinNegatCoord];
        end  % End of Width Threshold
    end % End of Looking for non empty values on both sides of the sulci
    Surfwidth.SurfData.faces(po,1:2) = [2*po-1 2*po];
end

%% ============ Refilling Sulcal Spam Map Neighoborhoods =============== %%

% Removing sulcal widths greater than a specified threshold
Surfwidth = Reorg_Surf(Surfwidth);
ind2del= find(sulcalWidth > opts.maxsulcalwidth);
if ~isempty(ind2del)
    sulcalWidth(ind2del) = 0;
end
if sum(sulcalWidth)
    Surfwidth.SurfData.vertices(2*ind2del-1,:) = zeros(length(ind2del),3);
    Surfwidth.SurfData.vertices(2*ind2del,:) = zeros(length(ind2del),3);
    diffVecs = Surfwidth.SurfData.vertices(1:2:end,:) - Surfo.SurfData.vertices; % Diference vector between the width line extreme and the surface
    norma = normm(diffVecs);
    dist2Surfo = sign(sum(diffVecs.*Normals,2)).*sqrt(sum(diffVecs.^2,2));
    
    % Width lines that intercept the median surface
    indcross = find(sign(dot(Surfwidth.SurfData.vertices(1:2:end,:) - Surfo.SurfData.vertices,Surfwidth.SurfData.vertices(2:2:end,:) - Surfo.SurfData.vertices,2)) == -1) ;
    dist2Surfo(indcross) = 0;
    
    
    % Removing sulcal widths that are far from mean surface
    ind2del= find(abs(dist2Surfo) > (opts.maxsulcalwidth)/2);
    if ~isempty(ind2del)
        sulcalWidth(ind2del) = 0;
        Surfwidth.SurfData.vertices(2*ind2del-1,:) = zeros(length(ind2del),3);
        Surfwidth.SurfData.vertices(2*ind2del,:) = zeros(length(ind2del),3);
        dist2Surfo(ind2del) = 0;
    end
    
    indsts = find(sulcalWidth == 0);
    
    indzeros = find(sum(Surfwidth.SurfData.vertices,2)==0);
    indface = find(sum(ismember(Surfwidth.SurfData.faces,indzeros),2)==2);
    Surfwidth.SurfData.faces(indface,:) = [];
    
    [Trip] = Vert_Neibp(double(Surfo.SurfData.faces),size(Surfo.SurfData.vertices,1),size(Surfo.SurfData.faces,1));
    Temp = sum(Trip);
    Trip(:,Temp==0) = [];
    temp = Trip(:,3:end);
    indz = find(temp == 0);
    temp(indz) = 1;
    
    sulcalWidth(indsts) = 0;
    dist2Surfo(indsts) = 0;
    indstsold = 0;
    while ~isempty(indsts)&sum(ismember(indstsold,indsts))~=length(indstsold);
        indstsold = indsts;
        %     iter = iter + 1;
        temp1 = sulcalWidth(temp);
        temp1(indz) = 0;
        temp1 = sum(temp1(indsts,:),2);
        
        ind2process = find(temp1 ~=0);
        %     ind2process = find(repmat(curvmap(Trip(indsts,2)),[1 size(temp1,2)]) - temp1 ~=0);
        
        % Recomputing Width Map
        Num = sulcalWidth(temp);
        Num(indz) = 0;
        Num = Num(indsts(ind2process),:);
        Den = sum(logical(Num),2);
        Num = sum(Num,2);
        sulcalWidth(indsts(ind2process)) = Num./(Den+eps); % Recomputing Width Map
        ind2del= find(sulcalWidth > opts.maxsulcalwidth);
        sulcalWidth(ind2del) = 0;
        
        % % %     % Recomputing Extreme Dist
        % % %     Num = dist2Surfo(temp);
        % % %     Num(indz) = 0;
        % % %     Num = Num(indsts(ind2process),:);
        % % %     Den = sum(logical(Num),2);
        % % %     Num = sum(Num,2);
        % % %     dist2Surfo(indsts(ind2process)) = Num./(Den+eps); % Recomputing Width Map
        % % %
        % % %     tempDistSurfo = Num./(Den+eps); % Recomputing Width Map
        % % %     Surfwidth.SurfData.vertices(2*indsts(ind2process)-1,:) = Surfo.SurfData.vertices(indsts(ind2process),:) + Normals(indsts(ind2process),:).*[tempDistSurfo tempDistSurfo tempDistSurfo];
        % % %
        % % %     tempSwidth = sulcalWidth(indsts(ind2process));
        % % %     Surfwidth.SurfData.vertices(2*indsts(ind2process),:) = Surfwidth.SurfData.vertices(2*indsts(ind2process)-1,:) + Normals(indsts(ind2process),:).*[tempSwidth tempSwidth tempSwidth];
        
        indsts(ind2process) = [];
        indsts = find(sulcalWidth == 0);
    end
end
Surfo.Is = sulcalWidth(:,1);
%% ========== End of Refilling Sulcal Spam Map Neighoborhoods ========== %%

%% ===================== End of Computing Sulcal Spam =================== %
% Outputs
varargout{1} = Surfwidth;
varargout{2} = sulcalWidth;
if opts.verbose
    toc;
end
return;