function varargout = Compute_Sulci_Metrics_V1(varargin);
%
% Syntax :
%  [Sulcmetrics, Surfo, SurfL] = Compute_Sulci_Metrics(SulcSurfMat, PialSurfMat, HullSurfMat, opts);
%
% This script computes the intercepting curve between two surfaces.
%
% Input Parameters:
%       SulcSurfMat             : Sulci Surface in matlab format
%       PialSurfMat             : Hemispheric Pial Surface in Matlab Format
%       HullSurfMat             : Hemispheric Hull Surface in Matlab Format
%       opts                    : Options
%                              opts.maxsulcawidth (mm). Maximum Sulcal Width
%                              opts.mindist2hull (mm). Minimum distance to the hull surface mm
%                              opts.angthreshold (degrees). Minimum angle allowed between normals
%
%
%
% Output Parameters:
%       SulcMetrics             : Sulcus metrics
%        Surfo                  : Reparametrized Sulcus
%        SurfL                  : Curves: Topline, Bottom line, Length,
%                                 Depth and width curves
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin < 4
    opts.maxsulcawidth = 12; % mm. Maximum Sulcal Width
    opts.mindist2hull = 5; % mm. Minimum distance to the hull surface mm
    opts.angthreshold = 5; % degrees. Minimum angle allowed between normals
    opts.ncurvp = 40; % Number of curve points
else
    opts = varargin{4};
end

if nargin < 3
    error('Three Inputs are needed');
    return
end
SulcSurfMat = varargin{1};
PialSurfMat = varargin{2};
HullSurfMat = varargin{3};
%% ==================== End of Input parameters  =========================%
% % % % % load('/media/COSAS/scripts/matlab.mat');

%% =========================== Main Program ============================= %

% Computing Maximum Edge length
maxdist =             sqrt(max(sum((SulcSurfMat.SurfData.vertices(SulcSurfMat.SurfData.faces(:,1),:)-SulcSurfMat.SurfData.vertices(SulcSurfMat.SurfData.faces(:,2),:)).^2,2)));
maxdist = max(maxdist,sqrt(max(sum((SulcSurfMat.SurfData.vertices(SulcSurfMat.SurfData.faces(:,1),:)-SulcSurfMat.SurfData.vertices(SulcSurfMat.SurfData.faces(:,3),:)).^2,2))));
maxdist = max(maxdist,sqrt(max(sum((SulcSurfMat.SurfData.vertices(SulcSurfMat.SurfData.faces(:,2),:)-SulcSurfMat.SurfData.vertices(SulcSurfMat.SurfData.faces(:,3),:)).^2,2))));


% Hull Surface planes
VertH = HullSurfMat.SurfData.vertices;FacesH = HullSurfMat.SurfData.faces;NormalsH = HullSurfMat.SurfData.VertexNormals;
nh = cross(VertH(FacesH(:,1),:)-VertH(FacesH(:,2),:),VertH(FacesH(:,3),:)-VertH(FacesH(:,2),:),2); DH = -1*dot(nh,VertH(FacesH(:,1),:),2);

%% ==================== Sulcus Walls Labelling ========================== %
SulcSurfMat = Sulcal_Face_Labelling(SulcSurfMat);
Surf2Process = SulcSurfMat;

%% ================== End of Sulcus Walls Labelling ===================== %


% ----------------------------------------------------------- %
%% =============== Computing Sulcal Length ============================= %%

disp(' ');
disp('Extracting the interception curve between the WM Hull Surface and Sulcus median mesh ...');
tic;
[sulchullIntercep, edgeLabels] = Intercept_Surface_with_Surface(HullSurfMat,  Surf2Process);


topCoords = Compute_mean_between_walls(sulchullIntercep, edgeLabels);


toc;
if size(sulchullIntercep,1) < 5 % At least 5 interception points are needed to construct a correct interception curve
    sulchullIntercep = '';
end
if ~isempty (sulchullIntercep)
%     topCoords = fitCurveTo3DPts(sulchullIntercep(:,1:3), sulchullIntercep(maxdistPair(1),1:3), sulchullIntercep(maxdistPair(2),1:3),size(sulchullIntercep,1)*2, 0 ); % Fitting a 3D curve
    
    
    disp(' ');
    disp('Creating perpendicular planes to the interception curve... ');
    tic;
    
    % ------- Creating perpendicular planes to the interception line (sulchullIntercep(:,1:3) or topCoords)
    P1 = [topCoords(1:size(topCoords,1)-2,1:3)];
    P2 = [topCoords(2:size(topCoords,1)-1,1:3)];
    P3 = [topCoords(3:size(topCoords,1),1:3)];
    n = cross(P1-P2,P3-P2,2); D = -1*dot(n,P1,2); %  Planes equations for each 3 points of the line
    t = P2-P1; u = P3-P1; v = P3-P2;
    t2 = sum((t.^2)')'; u2 = sum((u.^2)')'; w2 = sum((n.^2)')';
    c = P1+(repmat(t2.*dot(u,v,2),[1 size(u,2)]).*u-repmat(u2.*dot(t,v,2),[1 size(t,2)]).*t)./(2*repmat(w2,[1 size(P1,2)])+eps); % Center of the circunference
    r = 1/2*sqrt(t2.*u2.*dot(v,v,2)./w2); % Radius of the circunference
    
    % ------- Removing points with high curvature
    curv = ones(length(r),1)./r;
    ind = find(curv > 2);
    ind = ind +1;
    topCoords(ind,:) = [];
    iter = 0;
    while ~isempty(ind)| iter ==3
        iter = iter + 1;
        P1 = [topCoords(1:size(topCoords,1)-2,1:3)];
        P2 = [topCoords(2:size(topCoords,1)-1,1:3)];
        P3 = [topCoords(3:size(topCoords,1),1:3)];
        n = cross(P1-P2,P3-P2,2); D = -1*dot(n,P1,2);
        t = P2-P1; u = P3-P1; v = P3-P2;
        t2 = sum((t.^2)')'; u2 = sum((u.^2)')'; w2 = sum((n.^2)')';
        c = P1+(repmat(t2.*dot(u,v,2),[1 size(u,2)]).*u-repmat(u2.*dot(t,v,2),[1 size(t,2)]).*t)./(2*repmat(w2,[1 size(P1,2)])+eps);
        r = 1/2*sqrt(t2.*u2.*dot(v,v,2)./w2+eps);
        curv = ones(length(r),1)./r;
        ind = find(curv > 2);
        ind = ind +1;
        topCoords(ind,:) = [];
    end
    %topCoords = fitCurveTo3DPts(topCoords(:,1:3), topCoords(1,1:3), topCoords(end,1:3),size(topCoords,1)+2, 0 ); % Fitting a 3D curve
    
    
    % ------- Creating perpendicular planes to the interception line
    P1 = [topCoords(1:size(topCoords,1)-2,1:3)];
    P2 = [topCoords(2:size(topCoords,1)-1,1:3)];
    P3 = [topCoords(3:size(topCoords,1),1:3)];
    n = cross(P1-P2,P3-P2,2); D = -1*dot(n,P1,2); %  Planes equations for each 3 points of the line
    t = P2-P1; u = P3-P1; v = P3-P2;
    t2 = sum((t.^2)')'; u2 = sum((u.^2)')'; w2 = sum((n.^2)')';
    c = P1+(repmat(t2.*dot(u,v,2),[1 size(u,2)]).*u-repmat(u2.*dot(t,v,2),[1 size(t,2)]).*t)./(2*repmat(w2,[1 size(P1,2)])+eps); % Center of the circunference
    r = 1/2*sqrt(t2.*u2.*dot(v,v,2)./w2); % Radius of the circunference
    n = cross(P1-c,P2-c,2); D = -1*dot(n,c,2); % Perpendicular panes of the circunference
    n1 = cross(P2-c,n,2); D1 = -1*dot(n1,c,2); %n1 = [n1(1,:);n1;n1(end,:)];D1 = [D1(1);D1;D1(end)];
    topCoords = topCoords(2:end-1,:);
    toc;
    
    %% =================== Computing Top and Bottom Lines =============== %
    disp(' ');
    disp('Computing Top and Bottom lines ... ');
    cont = 0;
    tic;
    
    % First depth estimation
    Np = size(n1,1); % Number of perpendicular planes
    cont2del = 0;
    Plot_Surf(SulcSurfMat);
    for planes = 1:Np
        [intercLine, edgeLabels] = Intercept_Plane_with_Surface([n1(planes,:) D1(planes,:)],Surf2Process);
        

     
%         line(Xline,Yline,Zline,'Linewidth',2,'Color',[1 1 1]);
       
        if ~isempty(intercLine)
            
            topCoordsC = Compute_mean_between_walls(intercLine, edgeLabels);
            dist2Points = sqrt(sum((topCoordsC(:,1:3) - repmat(topCoords(planes,1:3),[size(topCoordsC,1) 1])).^2,2));
            [~,loc] = min(dist2Points);
            clustId = topCoordsC(loc,4);
            ind2del = find(topCoordsC(:,4) ~= clustId); % Indexes of the interception line that will be deleted
            topCoordsC(ind2del,:) = [];
            
            if ~isempty(topCoordsC)
%                 line(Xline,Yline,Zline,'Linewidth',2,'Color',[1 1 1]);
                cont = cont + 1;
                tempExtrem = [topCoordsC(1,1:3) ; topCoordsC(end,1:3)];
                
                %% ========== Detecting Both Sides of the Line ========== %
                intop = inpolyhedron(HullSurfMat.SurfData, tempExtrem);     % Verifying if all top lines are outside the hull
                %varargout = Plot_Plane([n1(planes,:) D1(planes,:)],topCoords(planes,:),3);
                caseopt = sum(intop);
                if caseopt == 1 % One point at each side of the hull
                    botpoint  = tempExtrem(intop == 1,:);
                    toppoint  = tempExtrem(intop == 0,:);
                elseif caseopt == 0 % Both points are outside the hull
                    tempc = [topCoords(planes,:);tempExtrem];
                    tempd = dist(tempc');
                    [ord,pos] = sort(tempd(1,2:end));
                    
                    botpoint  = tempExtrem(pos(1),:);        % Bottom point
                    toppoint  = tempExtrem(pos(2),:);   % Top point
                elseif  caseopt == 2 % Both points are inside the hull
                    tempc = [topCoords(planes,:);tempExtrem];
                    tempd = dist(tempc');
                    [ord,pos] = sort(tempd(1,2:end));
                    
                    botpoint  = tempExtrem(pos(2),:);        % Bottom point
                    toppoint  = tempExtrem(pos(1),:);   % Top point
                end
                hold on;
                plot3(topCoordsC(:,1),topCoordsC(:,2),topCoordsC(:,3),'-w','Linewidth',2,'Markersize',20);
                plot3(botpoint(:,1),botpoint(:,2),botpoint(:,3),'.y','Markersize',40); 
                plot3(toppoint(:,1),toppoint(:,2),toppoint(:,3),'.c','Markersize',40)
                %% ========== Detecting Both Sides of the Line ========== %
                % Sorting Boundary Lines
                if cont > 2
                    
                    % Top Line
                    intPoints = Intercept_Plane_with_LinesSegment([n1(planes,:) D1(planes,:)], topline(1:end-1,:),topline(2:end,:));
                    if ~isempty(intPoints)
                        
                        Nint = size(intPoints,1);
                        temp = zeros(size(topline,1)+Nint ,3);
                        [~,ord] = sort(intPoints(:,4));
                        intPoints = intPoints(ord,:);
                        temp(intPoints(:,4)+ [1:Nint]',:) = intPoints(:,1:3);
                        temp(setdiff(1:size(temp,1), intPoints(:,4)+ [1:Nint]')',:) = topline ;
                        topline = temp;
                        
                        %topline = [topline(1:intPoints(4),:) ; toppoint; topline(intPoints(4)+1:end,:)];
                    else
                        topline = [topline;toppoint]; % Top lines coordinates
                    end
                    intPoints = Intercept_Plane_with_LinesSegment([n1(planes,:) D1(planes,:)], bottomline(1:end-1,:),bottomline(2:end,:));
                    if ~isempty(intPoints)
                        
                        Nint = size(intPoints,1);
                        temp = zeros(size(bottomline,1)+Nint ,3);
                        [~,ord] = sort(intPoints(:,4));
                        intPoints = intPoints(ord,:);
                        temp(intPoints(:,4)+ [1:Nint]',:) = intPoints(:,1:3);
                        temp(setdiff(1:size(temp,1), intPoints(:,4)+ [1:Nint]')',:) = bottomline ;
                        bottomline = temp;
                        
                        %bottomline = [bottomline(1:intPoints(4),:) ; botpoint; bottomline(intPoints(4)+1:end,:)];
                    else
                        bottomline = [bottomline;botpoint]; % Top lines coordinates
                    end
                    
                elseif cont == 2
                    
                    % Top Line
                    topline = [topline;toppoint]; % Top lines coordinates
                    if sign(dot(topline(2,:) - topline(1,:), topCoordsC(2,1:3) - topCoordsC(1,1:3),2)) == -1;
                        topline = flipdim(topline,1);
                    end
                    
                    % Bottom Line
                    bottomline = [bottomline;botpoint]; % Top lines coordinates
                    if sign(dot(bottomline(2,:) - bottomline(1,:), topCoordsC(2,1:3) - topCoordsC(1,1:3),2)) == -1;
                        bottomline = flipdim(bottomline,1);
                    end
                    
                elseif cont == 1
                    topline    = toppoint; % Top lines coordinates
                    bottomline = botpoint; % Bottom lines coordinates
                end
            else
                cont2del =  cont2del + 1;
                delindex(cont2del) = planes;
            end
        end
        
        
    end
    if exist('delindex','var')
        topCoords(delindex,:) = [];
    end
    try
        tempSlength = max(max(dist(topline')));% Temporal Sulcal Length Estimation
        opts.ncurvp = 2*ceil((2/3)*tempSlength); % Number of points according to the sulcal length
        if opts.ncurvp <3
            opts.ncurvp = 5;
        end
        toplineCoords = fitCurveTo3DPts(topline, topline(1,:), topline(end,:),opts.ncurvp, 0 ); % Fitting a 3D curve
        bottomlineCoords = fitCurveTo3DPts(bottomline, bottomline(1,:), bottomline(end,:),opts.ncurvp, 0); % Fitting a 3D curve
        if size(toplineCoords,1) ~= size(bottomlineCoords,1)
            toplineCoords = Pts_fit_3Dline(toplineCoords, toplineCoords(1,:), toplineCoords(end,:),opts.ncurvp); % Fitting a 3D curve
            bottomlineCoords = Pts_fit_3Dline(bottomlineCoords, bottomlineCoords(1,:), bottomlineCoords(end,:),opts.ncurvp); % Fitting a 3D curve
        end
        
        intop = inpolyhedron(HullSurfMat.SurfData, toplineCoords);     % Verifying if all top lines are outside the hull
        indtop2rem = find(intop == 1);
        
        inbot = inpolyhedron(HullSurfMat.SurfData, bottomlineCoords);     % Verifying if all top lines are outside the hull
        indbot2rem = find(inbot == 0);
        ind2rem = unique([indtop2rem;indbot2rem]);
        if ~isempty(ind2rem)
            toplineCoords(ind2rem,:) = [];
            bottomlineCoords(ind2rem,:) = [];
            opts.ncurvp = size(bottomlineCoords,1);
        end
    catch
        toplineCoords = topline;
        bottomlineCoords = bottomline;
    end
    toc;
    %% ============ End of Computing Top and Bottom Lines =========== %
    
    %% =================== Computing Depth Lines ==================== %
    disp(' ');
    disp('Computing Depth lines ... ');
    cont = 0;
    tic;
    diffVectop = (toplineCoords(2:end,:) - toplineCoords(1:end-1,:)); % top line unitary vectors
    diffVectop = [diffVectop;diffVectop(end,:)];
    diffVectopbot = (bottomlineCoords(1:end,:) - toplineCoords(1:end,:)); % top-bottom unitary vectors
    crossVect = cross(diffVectop,diffVectopbot,2); % Vectorial product
    a = normm(crossVect);crossVect = crossVect./[a a a ];
    P1 = toplineCoords; % 1st point for Perpendicular plane creation
    P2 = bottomlineCoords; % 2nd point for Perpendicular plane creation
    P3 = toplineCoords + crossVect*10; % 3rd point for Perpendicular plane creation
    n1 = cross(P1-P2,P3-P2,2); D1 = -1*dot(n1,P2,2); %  Planes equations for each 3 points
    
    Nv = length(Surf2Process.SurfData.VertexNormals);
    
    SulcVert = [ 0 0 0 ];
    DlinesP = [0 0];
    DlinesF = [0 0];
    if opts.ncurvp > 10
        ndepth = round(opts.ncurvp/3); % Number of depth lines
    else
        ndepth = 8;
    end
    allobool=isnan(sum(n1,2));
    indinterp = find(allobool);
    indbase = find(allobool == 0);
    if ~isempty(indinterp)
        n1(indinterp,1) = interp1(indbase,n1(indbase,1),indinterp,'spline','extrap');
        n1(indinterp,2) = interp1(indbase,n1(indbase,2),indinterp,'spline','extrap');
        n1(indinterp,3) = interp1(indbase,n1(indbase,3),indinterp,'spline','extrap');
    end
    
    
    for planes = 1:size(n1,1)
        [intercLine, edgeLabels] = Intercept_Plane_with_Surface([n1(planes,:) D1(planes,:)],Surf2Process);
        topCoordsC = Compute_mean_between_walls(intercLine, edgeLabels);

        DepthCoords = fitCurveTo3DPts([toplineCoords(planes,:);topCoordsC(:,1:3);bottomlineCoords(planes,:)], toplineCoords(planes,:), bottomlineCoords(planes,:),ndepth, 0); % Fitting a 3D curve
        
        if size(DepthCoords,1) ~= ndepth
            DepthCoords = Pts_fit_3Dline([toplineCoords(planes,:);topCoordsC(:,1:3);bottomlineCoords(planes,:)], toplineCoords(planes,:), bottomlineCoords(planes,:),ndepth); % Fitting a 3D curve
        end
        SulcMetrics.depth.Coords{planes} = DepthCoords; % Depth Coordinates
        SulcVert = [SulcVert;DepthCoords];
        if (planes == 1)|(planes == opts.ncurvp)
            DlinesF = [DlinesF;[(1:(ndepth-1))' (2:ndepth)'] + ndepth*(planes-1)];
        else
            DlinesP = [DlinesP;[(1:(ndepth-1))' (2:ndepth)'] + ndepth*(planes-1)];
        end
        DepthProfile(planes,1) = sum(sqrt(sum(((DepthCoords(1:end-1,:) - DepthCoords(2:end,:)).^2)')));
    end
    for planes = 1:ndepth
        Temp = [planes:ndepth:(ndepth*(opts.ncurvp-1)+planes)]';
        if (planes >2) | (planes < ndepth)
            DlinesP = [DlinesP;[Temp(1:end-1) Temp(2:end)]];
        end
        templine = SulcVert(Temp(2:end),:);
        LengthProfile(planes,1) = sum(sqrt(sum(((templine(1:end-1,:) - templine(2:end,:)).^2)')));
    end
    SulcVert(1,:) = [];
    DlinesP(1,:) = [];
    DlinesF(1,:) = [];
    
    % Sulcal Length
    
    % Computing Interception Length
    verts = sulchullIntercep(:,1:3);
    intercpLength = sqrt(sum((verts(1:end-1,:) - verts(2:end,:)).^2,2));
    
    % Metrics
    SulcMetrics.length.measures = [sum(intercpLength); max(LengthProfile); min(LengthProfile); mode(LengthProfile) ; median(LengthProfile) ;  mean(LengthProfile);  std(LengthProfile)];
    SulcMetrics.length.profile = LengthProfile;
    
    % Sulcal Depth
    SulcMetrics.depth.measures = [max(DepthProfile); min(DepthProfile); mode(DepthProfile) ; median(DepthProfile) ;  mean(DepthProfile);  std(DepthProfile)];
    SulcMetrics.depth.profile = DepthProfile(:);
    
    toc;
    %% ================= Creating a New Sulci Surface (I) =========== %
    
    % Top Surface
    Surfo.SurfData.vertices = SulcVert;
    
    % First face indexes
    temp = repmat([2:ndepth-1],[2 1]);
    fac1 = [1;temp(:);ndepth];
    
    % Second face indexes
    fac2 = repmat([ndepth+1:2*ndepth-1],[2 1]);fac2= fac2(:);
    
    % Third face indexes
    fac3 = fac1*0;
    fac3(1:2:end) = fac1(1:2:end)+1;
    fac3(2:2:end) = fac2(2:2:end)+1;
    
    temp = repmat(ndepth*[1:opts.ncurvp-2],[length(fac1) 1]);
    Surfo.SurfData.faces = repmat([fac1 fac2 fac3],[opts.ncurvp-1 1]) + repmat([zeros(length(fac1),1);temp(:)],[1 3]);
    
    % Computing Normals
    Surfo = Compute_Surface_Normals(Surfo);
    Normals = Surfo.SurfData.VertexNormals;
    
    % Creating New Interception Line
    % [sulchullIntercep] = Intercept_Surface_with_Surface(HullSurfMat,  Surfo);
    
    %% ========= End of Creating a New Sulci Surface (I) =========== %
    
    %% ================= Computing Sulcal Spam  ===================== %
    disp(' ');
    disp('Computing Width lines ... ');
    cont = 0;
    tic;
    % Pial Surface planes
    VertP = PialSurfMat.SurfData.vertices;FacesP = PialSurfMat.SurfData.faces;NormalsP = PialSurfMat.SurfData.VertexNormals;
    np = cross(VertP(FacesP(:,1),:)-VertP(FacesP(:,2),:),VertP(FacesP(:,3),:)-VertP(FacesP(:,2),:),2); Dp = -1*dot(np,VertP(FacesP(:,1),:),2);
    
    P1 = Surfo.SurfData.vertices+opts.maxsulcawidth*Normals; P2 = Surfo.SurfData.vertices-opts.maxsulcawidth*Normals; % Creating perpendicular lines to sulcus walls
    Npoints = length(P1);
    cont = 0;
    Verts = Surfo.SurfData.vertices;
    
    if Npoints >200
        P = round(linspace(1,Npoints,200));
        P1 = P1(P,:);
        P2 = P2(P,:);
        Verts = Verts(P,:);
    end
    Npoints = length(P1);
    for po =1:Npoints
        Num = dot(np,repmat(P1(po,:),[size(np,1) 1]),2)+Dp; Den = dot(np,repmat(P2(po,:)-P1(po,:),[size(np,1) 1]),2); % Creating Lines
        
        t = -1*Num./Den;clear Num Den; % Lines parameters
        xint = single(repmat(P1(po,1)',[size(FacesP,1) 1])+t.*(repmat((P2(po,1)-P1(po,1))',[size(FacesP,1) 1]))); % Line parametrization
        yint = single(repmat(P1(po,2)',[size(FacesP,1) 1])+t.*(repmat((P2(po,2)-P1(po,2))',[size(FacesP,1) 1])));
        zint = single(repmat(P1(po,3)',[size(FacesP,1) 1])+t.*(repmat((P2(po,3)-P1(po,3))',[size(FacesP,1) 1])));clear t;
        
        PpP1 =  VertP(FacesP(:,1),:)-[xint yint zint];
        PpP2 =  VertP(FacesP(:,2),:)-[xint yint zint];
        PpP3 =  VertP(FacesP(:,3),:)-[xint yint zint];
        angP2p3 = (acos(dot(PpP2,PpP3,2)./(normm(PpP2).*normm(PpP3))))*180/pi; % Angles between each face point and the interception point inpial surface
        angP1p3 = (acos(dot(PpP1,PpP3,2)./(normm(PpP1).*normm(PpP3))))*180/pi;
        angP1p2 = (acos(dot(PpP1,PpP2,2)./(normm(PpP1).*normm(PpP2))))*180/pi;
        
        ind = find(round(angP2p3+angP1p3+angP1p2) == 360);
        [inte,indord] = unique([xint(ind) yint(ind) zint(ind)],'rows');ind = ind(indord); % Line interception with pial surface
        
        
        t = NormalsP(FacesP(ind,:),:);
        intNormals = reshape(mean(reshape(reshape(t,[length(ind) prod(size(t))/length(ind)])',[3 length(ind)*3]))',[3 length(ind)])';
        norma = normm(intNormals);
        intNormals = intNormals./[norma norma norma];
        
        
        if ~isempty(inte)
            %     orientsign = sign((dot(inte -repmat(Surf2Process.SurfData.vertices(po,:),[size(inte,1) 1]),repmat(Surf2Process.SurfData.VertexNormals(po,:),[size(inte,1) 1]),2))); % Detecting Orientation
            distances = sqrt(sum((inte - repmat(Verts(po,:),[size(inte,1) 1])).^2,2)); % Computing distance to the pial surface
            
            [valMin,posit] = min(distances);
            if valMin < opts.maxsulcawidth
                MinPositCoord = inte(posit,:);
                
                
                t1 = inte - repmat(MinPositCoord,[size(inte,1) 1]);
                t2 = repmat(intNormals(posit,:),[size(inte,1) 1]);
                
                indPositDir = find(sum(t1.*t2,2) > 0);
                if ~isempty(indPositDir)
                    distances = sqrt(sum((inte(indPositDir,:) - repmat(MinPositCoord,[length(indPositDir) 1])).^2,2));
                    ind = find((distances == min(nonzeros(distances))));
                    MinNegatCoord = inte(indPositDir(ind),:);
                    sulcw = distances(ind);
                    
                    %                 distances = sqrt(sum((inte - repmat(MinPositCoord,[size(inte,1) 1])).^2,2));
                    %                 ind = find((distances == min(nonzeros(distances))));
                    %                 MinNegatCoord = inte(ind,:);
                    %                 sulcw = distances(ind);
                    if sulcw < opts.maxsulcawidth  % Width Threshold
                        cont = cont  + 1;
                        sulcalwidth(cont,:) = [po sulcw];
                        %plot3([MinPositCoord(1,1);MinNegatCoord(1,1)],[MinPositCoord(1,2);MinNegatCoord(1,2)],[MinPositCoord(1,3);MinNegatCoord(1,3)],'-','Color',[0 0 1],'Markersize',30);
                        SurfL(5).SurfData.vertices(2*cont-1:2*cont,:) = [MinPositCoord;MinNegatCoord];
                        SurfL(5).SurfData.faces(cont,1:2) = [2*cont-1 2*cont];
                    end  % End of Width Threshold
                end % End of Looking for non empty values on both sides of the sulci
            end
        end
    end
    %
    SulcMetrics.width.measures = [max(sulcalwidth(:,2)); min(sulcalwidth(:,2)) ; mode(sulcalwidth(:,2)) ; median(sulcalwidth(:,2))  ;  mean(sulcalwidth(:,2));  std(sulcalwidth(:,2))];
    SulcMetrics.width.profile = sulcalwidth(:,2);
    toc;
    %% ================= End of Computing Sulcal Spam ============== %
    
    
    
    
    %% == Detecting Sign to create the parametrized sulcus surface == %
    cont = 0;
    % Pial Surface planes
    VertS = SulcSurfMat.SurfData.vertices;FacesS = SulcSurfMat.SurfData.faces;
    np = cross(VertS(FacesS(:,1),:)-VertS(FacesS(:,2),:),VertS(FacesS(:,3),:)-VertS(FacesS(:,2),:),2); Dp = -1*dot(np,VertS(FacesS(:,1),:),2);
    
    P1 = Surfo.SurfData.vertices+100*Normals; P2 = Surfo.SurfData.vertices-100*Normals; % Creating perpendicular lines to sulcus walls
    Npoints = length(P1);
    cont = 0;
    Verts = Surfo.SurfData.vertices;
    
    if Npoints >200
        P = round(linspace(1,Npoints,200));
        P1 = P1(P,:);
        P2 = P2(P,:);
        Verts = Verts(P,:);
    end
    Npoints = length(P1);
    for po =1:Npoints
        Num = dot(np,repmat(P1(po,:),[size(np,1) 1]),2)+Dp; Den = dot(np,repmat(P2(po,:)-P1(po,:),[size(np,1) 1]),2); % Creating Lines
        
        t = -1*Num./Den;clear Num Den; % Lines parameters
        xint = single(repmat(P1(po,1)',[size(FacesS,1) 1])+t.*(repmat((P2(po,1)-P1(po,1))',[size(FacesS,1) 1]))); % Line parametrization
        yint = single(repmat(P1(po,2)',[size(FacesS,1) 1])+t.*(repmat((P2(po,2)-P1(po,2))',[size(FacesS,1) 1])));
        zint = single(repmat(P1(po,3)',[size(FacesS,1) 1])+t.*(repmat((P2(po,3)-P1(po,3))',[size(FacesS,1) 1])));clear t;
        
        PpP1 =  VertS(FacesS(:,1),:)-[xint yint zint];
        PpP2 =  VertS(FacesS(:,2),:)-[xint yint zint];
        PpP3 =  VertS(FacesS(:,3),:)-[xint yint zint];
        angP2p3 = (acos(dot(PpP2,PpP3,2)./(normm(PpP2).*normm(PpP3))))*180/pi; % Angles between each face point and the interception point inpial surface
        angP1p3 = (acos(dot(PpP1,PpP3,2)./(normm(PpP1).*normm(PpP3))))*180/pi;
        angP1p2 = (acos(dot(PpP1,PpP2,2)./(normm(PpP1).*normm(PpP2))))*180/pi;
        
        ind = find(round(angP2p3+angP1p3+angP1p2) == 360);
        [inte,indord] = unique([xint(ind) yint(ind) zint(ind)],'rows');ind = ind(indord); % Line interception with pial surface
        if ~isempty(inte)
            orientsign = sign((dot(inte -repmat(Surfo.SurfData.vertices(po,:),[size(inte,1) 1]),repmat(Normals(po,:),[size(inte,1) 1]),2))); % Detecting Orientation
            distances = sqrt(sum((inte - repmat(Verts(po,:),[size(inte,1) 1])).^2,2)); % Computing distance to the pial surface
            [maxValDist(po),loc] = max(distances); % Detecting Median mesh width  
             meanOrient(po)= orientsign(loc); %orientation to create the parametrized surface
        end
    end
    %% == Detecting Sign to create the parametrized sulcus surface == %
    
    
    
    
    %% ================= Creating a New Sulci Surface (II) =========== %
    % Bottom Surface
    clear Surft;
    Surfo.SurfData.vertices = [Surfo.SurfData.vertices;Surfo.SurfData.vertices + mode(meanOrient)*Normals*mean(maxValDist/2)];
    Surfo.SurfData.faces = [Surfo.SurfData.faces;flipdim(Surfo.SurfData.faces,2)+max(Surfo.SurfData.faces(:))];
    Surfo.SurfData.VertexNormals = [Surfo.SurfData.VertexNormals;-1*Surfo.SurfData.VertexNormals];
    
    indexes = [1:ndepth:(ndepth*(opts.ncurvp -1) +1) (ndepth*(opts.ncurvp -1) +2):(ndepth*opts.ncurvp-1)  ndepth*opts.ncurvp:-1*ndepth:ndepth (ndepth-1):-1:1]';
    
    toplineindexes = [1:ndepth:(ndepth*(opts.ncurvp -1) +1)];
    toplineindexes = [toplineindexes(1:end-1)' toplineindexes(2:end)'];
    bottomlineindexes = [ndepth:ndepth:ndepth*opts.ncurvp];
    bottomlineindexes = [bottomlineindexes(1:end-1)' bottomlineindexes(2:end)'];
    
    indexes = [indexes;indexes + ndepth*opts.ncurvp];
    
    nindexes = length(indexes)/2;
    % First face index
    temp = repmat([2:nindexes-1],[2 1]);
    fac1 = [1;temp(:);nindexes];
    
    % Second face index
    fac2 = repmat([nindexes+1:2*nindexes-1],[2 1]);fac2= fac2(:);
    
    % Third face index
    fac3 = fac1*0;
    fac3(1:2:end) = fac1(1:2:end)+1;
    fac3(2:2:end) = fac2(2:2:end)+1;
    nfaces = [indexes(fac1) indexes(fac2)  indexes(fac3) ] ;
    Surfo.SurfData.faces = [Surfo.SurfData.faces;nfaces];
    FV2=smoothpatch(Surfo.SurfData,1,2);
    Surfo.SurfData.vertices = FV2.vertices;clear FV2;
    
    SulcVert = Surfo.SurfData.vertices(1:length(Surfo.SurfData.vertices)/2,:);
    
    
    % Saving Line Surfaces
    % Saving Plane Lines
    temp = unique(DlinesP);
    neworder = [1:length(temp)]';
    [a,t] = ismember(DlinesP,temp);
    vert2delet = find(ismember([1:length(SulcVert)],unique(DlinesP(:))) == 0);
    Surft.SurfData.vertices = SulcVert;
    Surft.SurfData.vertices(vert2delet,:) = [];
    Surft.SurfData.faces = t;
    SurfL(1) = Surft;
    
    % Saving Lateral Boundaries
    temp = unique(DlinesF);
    neworder = [1:length(temp)]';
    [a,t] = ismember(DlinesF,temp);
    vert2delet = find(ismember([1:length(SulcVert)],unique(DlinesF(:))) == 0);
    Surft.SurfData.vertices = SulcVert;
    Surft.SurfData.vertices(vert2delet,:) = [];
    Surft.SurfData.faces = t;
    SurfL(2) = Surft;
    
    % Saving Top Line
    temp = unique(toplineindexes);
    neworder = [1:length(temp)]';
    [a,t] = ismember(toplineindexes,temp);
    vert2delet = find(ismember([1:length(SulcVert)],unique(toplineindexes(:))) == 0);
    Surft.SurfData.vertices = SulcVert;
    Surft.SurfData.vertices(vert2delet,:) = [];
    Surft.SurfData.faces = t;
    SurfL(3) = Surft;
    %
    % Saving Bottom Line
    temp = unique(bottomlineindexes);
    neworder = [1:length(temp)]';
    [a,t] = ismember(bottomlineindexes,temp);
    vert2delet = find(ismember([1:length(SulcVert)],unique(bottomlineindexes(:))) == 0);
    Surft.SurfData.vertices = SulcVert;
    Surft.SurfData.vertices(vert2delet,:) = [];
    Surft.SurfData.faces = t;
    SurfL(4) = Surft;
    
    % Saving Interception Line
    Surft.SurfData.vertices = sulchullIntercep(:,1:3);
    Nintvert = size(Surft.SurfData.vertices,1);
    Surft.SurfData.faces = [[1:Nintvert-1]' [2:Nintvert]'];
    SurfL(6) = Surft;
    
    %% =============== End of Creating a New Surface =============== %%
    %% ================= End of Computing Depth Lines =================== %
else
    % Reparametrized Surface
    Surfo = '';
    
    % Curves
    SurfL = '';
    
    % Sulcal Width
    SulcMetrics.width.measures = zeros(1,6);
    SulcMetrics.width.profile = 0;
    % Sulcal Length
    SulcMetrics.length.measures = zeros(1,7);
    SulcMetrics.length.profile = 0;
    
    % Sulcal Depth
    SulcMetrics.depth.measures = zeros(1,6);
    SulcMetrics.depth.profile = 0;
end

varargout{1} = SulcMetrics;
varargout{2} = Surfo;
varargout{3} = SurfL;
return

%% ======================= Internal Functions =========================== %

function [cPts] = Pts_fit_3Dline(XYZ, stPt, endPt, Npoints);
MAX_ITERATIONS = 5;
completedIterations = 0;

% Initialise control pts linearly between start/end anchors
cPts = interp1(0:1, [stPt; endPt], linspace(0,1,Npoints));

while completedIterations < MAX_ITERATIONS
    % Get the nearest-neighbour cntrl pt for each of the sample points
    sqDists = cellfun(@(x)sum(bsxfun(@minus, XYZ, x).^2,2), num2cell(cPts,2),'UniformOutput', false);
    [~, nnIdxs] = min(cat(2,sqDists{:}),[],2);
    % Keep the anchors, update inner cPts to the mean of their linked input pts
    for i = 2:Npoints-1
        cPts(i,:) = mean(XYZ(nnIdxs==i,:),1);
    end
    % Handle any cPts that didn't have linked pts so their mean became NaN
    goodIdxs = find(~isnan(cPts(:,1)));
    badIdxs = find(isnan(cPts(:,1)));
    cPts(badIdxs,:) = interp1(goodIdxs, cPts(goodIdxs,:), badIdxs);
    % Re-spread the control points out linearly
    cPtCumSumDists = cumsum([0; sqrt(sum(diff(cPts,1).^2,2))]);
    cPts = interp1(cPtCumSumDists, cPts, linspace(0, cPtCumSumDists(end), Npoints));
    completedIterations = completedIterations + 1;
end
return;

function norma = normm(M)
norma = sqrt(sum((M').^2))';
return

%% =====================  Internal functions ============================ %
function topCoords = Compute_mean_between_walls(sulchullIntercep, edgeLabels);
[maxdistPair] = Sort_3D_pointCloud(sulchullIntercep);


Nc = size(maxdistPair,1);
topCoords = [0 0 0 0];
for i = 1:Nc
    ind = find(sulchullIntercep(:,4) == i);
    sulchullIntercepCl = sulchullIntercep(ind,:);
    edgeLabelsCl = edgeLabels(ind,:);
    tempLabelsCl = unique(edgeLabelsCl(:,1));
    maxdistPairCl = maxdistPair(i,:);
    for k = 1:length(tempLabelsCl)
        
        indpoints = find(edgeLabelsCl(:,1) == tempLabelsCl(k)&edgeLabelsCl(:,2) == tempLabelsCl(k));
        if ~isempty (sulchullIntercepCl)
            tempPoints = [sulchullIntercepCl(maxdistPairCl(1),1:3);sulchullIntercepCl(indpoints,1:3) ;sulchullIntercepCl(maxdistPairCl(2),1:3)];
            if k == 1
                topCoordsCl = fitCurveTo3DPts(tempPoints, tempPoints(1,:), tempPoints(end,:),size(sulchullIntercepCl,1)*2, 0 ); % Fitting a 3D curve
            else
                topCoordsCl = topCoordsCl + fitCurveTo3DPts(tempPoints, tempPoints(1,:), tempPoints(end,:),size(sulchullIntercepCl,1)*2, 0 ); % Fitting a 3D curve
            end
        end
        
    end
    topCoords = [topCoords;topCoordsCl/2 i*ones(size(topCoordsCl,1),1)];
end
topCoords(1,:) = [];
return