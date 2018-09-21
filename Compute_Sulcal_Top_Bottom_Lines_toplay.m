function varargout =  Compute_Sulcal_Top_Bottom_Lines_toplay(varargin);
%
% Syntax :
%       [bottomlineCoords, toplineCoords] = Compute_Sulcal_Top_Bottom_Lines(HullSurfMat, Surf2Process, topCoords, intNormals, opts);
%
% This function computes top and bottom lines that ride over and below the
% sulcal median mesh.
%
%
% Input Parameters:
%       HullSurfMat             : Hull Surface (Matlab Format)
%       Surf2Process            : Sulcal Median Mesh (Matlab Format).
%       topCoords               : Interception curve between Hull and Sulcal
%                                 median surfaces.
%       intNormals              : Normal vectors to the interception curve.
%       opts                    : Options:
%                                    - opts.maxsulcalwidth (Maximum Sulcal
%                                    Width in mm).
%                                    - opts.verbose (Verbose boolean
%                                    variable).
%                                    - opts.epsilon (Tolerance to estimate
%                                    barycentric coordinates).
%
% Output Parameters:
%
%      bottomlineCoords         : Coordinates of bottom lines.
%      toplineCoords            : Coordinates of top lines.
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% February 13th 2015
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin<4 % the indispensable input arguments are not provided
    error('Four inputs are mandatory');
    return
else
    HullSurfMat  = varargin{1};
    Surf2Process = varargin{2};
    topCoords    = varargin{3};
    intNormals   = varargin{4};
    
    HullSurfMat = Surface_Checking(HullSurfMat);
    Surf2Process = Surface_Checking(Surf2Process);
end
%% =========================== Input parameters  =========================%
if nargin < 5
    opts.maxsulcalwidth  = 12;     % mm. Maximum Sulcal Width
    opts.verbose = 1;              % Verbose boolean variable
    opts.epsilon = 10^-5;          % Tolerance to estimate barycentric coordinates
else
    opts = varargin{5};
    if ~isfield(opts,'maxsulcalwidth')
        opts.maxsulcalwidth = 12;      % Maximum Sulcal Width
    end
    if ~isfield(opts,'verbose')
        opts.verbose = 1;              % Verbose boolean variable
    end
    if ~isfield(opts,'epsilon')
        opts.epsilon = 1;              % Tolerance to estimate barycentric coordinates
    end
end
if nargin > 5
    error('To many inputs');
    return;
end

if nargout > 2
    error('To many outputs');
    return;
end
%% ==================== End of Input parameters  ======================== %

% load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/testmetrics.mat');

% ------- Creating perpendicular planes to the interception line
[T,N,B,k,t] = frenet(topCoords(:,1),topCoords(:,2),topCoords(:,3));

n1 = cross(intNormals,T,2);
n1 = cross(n1,intNormals,2);
D1 = -1*dot(n1,topCoords(:,1:3),2);

%% =================== Computing Top and Bottom Lines =============== %
if opts.verbose
    disp(' ');
    disp('Computing Top and Bottom lines ... ');
    tic;
end
% First depth estimation
% h =  Plot_Surf(Surf2Process);h(1).FaceAlpha = .5;hold on

% Hull Surface planes (each triangle is considered as a plane)
VertH = HullSurfMat.SurfData.vertices;FacesH = HullSurfMat.SurfData.faces;NormalsH = HullSurfMat.SurfData.VertexNormals;
np = cross(VertH(FacesH(:,1),:)-VertH(FacesH(:,2),:),VertH(FacesH(:,3),:)-VertH(FacesH(:,2),:),2); Dp = -1*dot(np,VertH(FacesH(:,1),:),2);

v0 = VertH(FacesH(:,3),:) - VertH(FacesH(:,1),:);
v1 = VertH(FacesH(:,2),:) - VertH(FacesH(:,1),:);
dot00 = dot(v0, v0, 2);
dot01 = dot(v0, v1, 2);
dot11 = dot(v1, v1, 2);
invDenom = 1 ./ (dot00 .* dot11 - dot01 .* dot01);


bottomline = [ 0 0 0]; % Bottom lines coordinates
topline    = [ 0 0 0]; % Top lines coordinates
Np = size(n1,1); % Number of perpendicular planes
for planes = 1:Np
    [intercLine, edgeLabels] = Intercept_Plane_with_Surface([n1(planes,:) D1(planes,:)],Surf2Process); % Interception between each plane and sulcal surface
    if ~isempty(intercLine)
        topCoordsC = Compute_mean_between_walls(intercLine); % Mean line of the interception line
        
        if length(unique(topCoordsC(:,4))) > 1
            indMultclust = find(diff(topCoordsC(:,4)));
            if ~isempty(indMultclust)
                dist2Points = sqrt(sum((topCoordsC(:,1:3) - repmat(topCoords(planes,1:3),[size(topCoordsC,1) 1])).^2,2));
                [~,loc] = min(dist2Points);
                clustId = topCoordsC(loc,4);
                if clustId ~=1
                    [topCoordsC] = flipdim(topCoordsC,1);
                    indMultclust = length(topCoordsC) - indMultclust;
                end
                
                % Detecting parallel interceptions
                aVector = (topCoordsC(indMultclust(1),1:3) - topCoordsC(1,1:3));
                aVector = repmat(aVector,[length(indMultclust) 1]);
                bVector = [topCoordsC(indMultclust+1,1) - topCoordsC(indMultclust(1),1) topCoordsC(indMultclust+1,2) - topCoordsC(indMultclust(1),2) topCoordsC(indMultclust+1,3) - topCoordsC(indMultclust(1),3)];
                
                angles = real(acos(sum(aVector.*bVector,2)./(sum(aVector.*aVector,2).*sum(bVector.*bVector,2))))*180/pi;
                clusIds = topCoordsC(indMultclust+1,4);
                clusIds2del = clusIds(angles>=20);
                ind2del = find(ismember(topCoordsC(:,4),clusIds2del)); % Indexes of the interception line that will be deleted
                topCoordsC(ind2del,:) = [];
%                 plot3(topCoordsC(:,1),topCoordsC(:,2),topCoordsC(:,3),'-r','Linewidth',2,'Markersize',20);
                
            end
        end
        if ~isempty(topCoordsC)
            tempExtrem = [topCoordsC(1,1:3) ; topCoordsC(end,1:3)];
            
            %% ========== Detecting Both Sides of the Line ========== %
            
            % Top - Bottom Vector
            P1  = tempExtrem(1,:);
            P2  = tempExtrem(2,:);
            tempNorm = normm(P1 - P2);
            vectOrient = [P1 - P2]./tempNorm;
            
            % Interception with Hull Surface planes
            Num = dot(np,repmat(P1,[size(np,1) 1]),2)+Dp; Den = dot(np,repmat(P2-P1,[size(np,1) 1]),2); % Creating Lines
            t = -1*Num./Den;clear Num Den; % Lines parameters
            xint = single(repmat(P1(1,1)',[size(FacesH,1) 1])+t.*(repmat((P2(1,1)-P1(1,1))',[size(FacesH,1) 1]))); % Line parametrization
            yint = single(repmat(P1(1,2)',[size(FacesH,1) 1])+t.*(repmat((P2(1,2)-P1(1,2))',[size(FacesH,1) 1])));
            zint = single(repmat(P1(1,3)',[size(FacesH,1) 1])+t.*(repmat((P2(1,3)-P1(1,3))',[size(FacesH,1) 1])));clear t;
            
            v2 = [xint yint zint]- VertH(FacesH(:,1),:);
            
            % Useful dot products
            dot02 = dot(v0, v2, 2);
            dot12 = dot(v1, v2, 2);
            
            % Obtaining barycentric coordinates 
            u = (dot11 .* dot02 - dot01 .* dot12) .* invDenom;
            v = (dot00 .* dot12 - dot01 .* dot02) .* invDenom;
            out = u >= -opts.epsilon & v >= -opts.epsilon & (u+v-opts.epsilon) <= 1;
            ind = find(out);
            
            [inte,indord] = unique([xint(ind) yint(ind) zint(ind)],'rows');ind = ind(indord); % Line interception with hull surface
            if ~isempty(inte)
                t = NormalsH(FacesH(ind,:),:);
                intNormals = reshape(mean(reshape(reshape(t,[length(ind) prod(size(t))/length(ind)])',[3 length(ind)*3]))',[3 length(ind)])';
                norma = normm(intNormals);
                intNormals = intNormals./[norma norma norma];
% %                 %     orientsign = sign((dot(inte -repmat(Surf2Process.SurfData.vertices,[size(inte,1) 1]),repmat(Surf2Process.SurfData.VertexNormals,[size(inte,1) 1]),2))); % Detecting Orientation
                distancesP1 = sqrt(sum((inte - repmat(P1,[size(inte,1) 1])).^2,2)); % Computing distance to the pial surface
                distancesP2 = sqrt(sum((inte - repmat(P2,[size(inte,1) 1])).^2,2)); % Computing distance to the pial surface
                distancesIC = sqrt(sum((inte - repmat(topCoords(planes,1:3),[size(inte,1) 1])).^2,2)); % Computing distance to the hull surface
% %                 
                 [~,posit] = min(distancesP1+distancesP2+distancesIC);
                
%                 MinPositCoord = inte(posit,:);
                MinPositNorma = intNormals(posit,:);
                
                if sum(MinPositNorma.*vectOrient,2) >=0
                    botpoint  = P2;
                    toppoint  = P1;
                else
                    botpoint  = P1;
                    toppoint  = P2;
                end
                bottomline = [bottomline;botpoint]; % Bottom lines coordinates
                topline    = [topline;toppoint];    % Top lines coordinates
                
            end
            %% ========== End of Detecting Both Sides of the Line ========== %
        else
            cont2del =  cont2del + 1;
            delindex(cont2del) = planes;
        end
    end
    %     hold on;
    %     plot3(topCoordsC(:,1),topCoordsC(:,2),topCoordsC(:,3),'-w','Linewidth',2,'Markersize',20);
%          plot3(botpoint(:,1),botpoint(:,2),botpoint(:,3),'.y','Markersize',40);
    %     plot3(toppoint(:,1),toppoint(:,2),toppoint(:,3),'.c','Markersize',40)
    %     Plot_Plane([n1(planes,:) D1(planes,:)],topCoords(planes,:),20)
end
bottomline(1,:) = [];
topline(1,:) = [];

locbottom = dsearchn(Surf2Process.SurfData.vertices,bottomline); % Closest bottom points to the surface
loctop    = dsearchn(Surf2Process.SurfData.vertices,topline);    % Closest top points to the surface


Graph = surface2graph(Surf2Process);
D = graphallshortestpaths(Graph);
tempVar = D(loctop,loctop);
[Xtop,Ytop] = find(tempVar == max(tempVar(:)));

tempVar = D(locbottom,locbottom);
[Xbot,Ybot] = find(tempVar == max(tempVar(:)));

if exist('delindex','var')
    topCoords(delindex,:) = [];
end

% Removing recursively bottomline points with high curvature 
[~,~,~,topCurv,~] = frenet(bottomline(:,1),bottomline(:,2),bottomline(:,3));
ind = find(topCurv>2);
while ~isempty(ind)
    bottomline(ind,:) = [];
    [~,~,~,topCurv,~] = frenet(bottomline(:,1),bottomline(:,2),bottomline(:,3));
    ind = find(topCurv>2);

end

% Removing recursively topline points with high curvature 
[~,~,~,topCurv,~] = frenet(topline(:,1),topline(:,2),topline(:,3));
ind = find(topCurv>2);
while ~isempty(ind)
    topline(ind,:) = [];
    [~,~,~,topCurv,~] = frenet(topline(:,1),topline(:,2),topline(:,3));
    ind = find(topCurv>2);

end

try
    tempSlength = max(max(dist(topline')));% Temporal Sulcal Length Estimation
    opts.ncurvp = 2*ceil((2/3)*tempSlength); % Number of points according to the sulcal length
    if opts.ncurvp <3
        opts.ncurvp = 5;
    end
    toplineCoords    = fitCurveTo3DPts(topline, topline(Xtop(1),:), topline(Ytop(1),:),opts.ncurvp, 0 ); % Fitting a 3D curve
    bottomlineCoords = fitCurveTo3DPts(bottomline, bottomline(Xbot(1),:), bottomline(Ybot(1),:),opts.ncurvp, 0 ); % Fitting a 3D curve
    if size(toplineCoords,1) ~= size(bottomlineCoords,1)
        toplineCoords = Pts_fit_3Dline(toplineCoords, toplineCoords(1,:), toplineCoords(end,:),opts.ncurvp); % Fitting a 3D curve
        bottomlineCoords = Pts_fit_3Dline(bottomlineCoords, bottomlineCoords(1,:), bottomlineCoords(end,:),opts.ncurvp); % Fitting a 3D curve
    end
catch
    toplineCoords = topline;
    bottomlineCoords = bottomline;
    opts.ncurvp = size(bottomlineCoords,1);
    
end

t1 = sqrt(sum((bottomlineCoords(1,:) - toplineCoords(1,:)).^2)) + sqrt(sum((bottomlineCoords(end,:) - toplineCoords(end,:)).^2));
t2 = sqrt(sum((bottomlineCoords(1,:) - toplineCoords(end,:)).^2)) + sqrt(sum((bottomlineCoords(end,:) - toplineCoords(1,:)).^2));
if t2 < t1
    bottomlineCoords = flipdim(bottomlineCoords,1);
end
bottomlineCoords = [smooth(bottomlineCoords(:,1)) smooth(bottomlineCoords(:,2)) smooth(bottomlineCoords(:,3))];
toplineCoords    = [smooth(toplineCoords(:,1)) smooth(toplineCoords(:,2)) smooth(toplineCoords(:,3))];

Npoints = size(bottomlineCoords,1);

% Removing recursively bottomlineCoords points with high curvature 
[~,~,~,topCurv,~] = frenet(bottomlineCoords(:,1),bottomlineCoords(:,2),bottomlineCoords(:,3));
ind = find(topCurv>2);
while ~isempty(ind)
    bottomlineCoords(ind,:) = [];
    [~,~,~,topCurv,~] = frenet(bottomlineCoords(:,1),bottomlineCoords(:,2),bottomlineCoords(:,3));
    ind = find(topCurv>2);

end

% Smoothing bottomlineCoords and obtain a curve with a specified number of
% points
[bottomlineCoords] = fitCurveTo3DPts(bottomlineCoords, bottomlineCoords(1,:), bottomlineCoords(end,:), Npoints, 0);


% Removing recursively toplineCoords points with high curvature 
[~,~,~,topCurv,~] = frenet(toplineCoords(:,1),toplineCoords(:,2),toplineCoords(:,3));
ind = find(topCurv>2);
while ~isempty(ind)
    toplineCoords(ind,:) = [];
    [~,~,~,topCurv,~] = frenet(toplineCoords(:,1),toplineCoords(:,2),toplineCoords(:,3));
    ind = find(topCurv>2);

end

% Smoothing toplineCoords and obtain a curve with a specified number of
% points
[toplineCoords] = fitCurveTo3DPts(toplineCoords, toplineCoords(1,:), toplineCoords(end,:), Npoints, 0);

if opts.verbose
    toc;
end
%% ============ End of Computing Top and Bottom Lines =================== %

% Outputs
varargout{1} = bottomlineCoords;
varargout{2} = toplineCoords;

return;