function varargout = Compute_Node_Metrics_toplay(varargin);
%
% Syntax :
%  [Sulcmetrics, Surfo, SurfL] = Compute_Node_Metrics(SulcSurfMat, PialSurfMat, HullSurfMat, opts);
%
% This script computes the depth, length and width for each sulcal node.
%
% Input Parameters:
%       SulcSurfMat             : Sulci Surface in matlab format
%       PialSurfMat             : Hemispheric Pial Surface in Matlab Format
%       HullSurfMat             : Hemispheric Hull Surface in Matlab Format
%       opts                    : Options
%                                 - opts.maxsulcalwidth (in mm. Maximum
%                                 Sulcal Width)
%                                 - opts.maxsulcaldepth (in mm. Maximum
%                                 Sulcal Depth).
%                                 - opts.maxsulcallength (mm. Maximum
%                                 Sulcal Length).
%                                 - opts.ncurvp (Number of curve points).
%                                 - opts.verbose (Number of curve points).
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
    opts.maxsulcalwidth  = 12; % mm. Maximum Sulcal Width
    opts.maxsulcaldepth  = 70; % mm. Maximum Sulcal Depth
    opts.maxsulcallength = 200; % mm. Maximum Sulcal Length
    opts.ncurvp  = 40; % Number of curve points
    opts.verbose = 0; % Number of curve points
    
else
    opts = varargin{4};
    if ~isfield(opts,'maxsulcalwidth')
        opts.maxsulcalwidth = 12;      % Maximum Sulcal Width
    end
    if ~isfield(opts,'maxsulcaldepth')
        opts.maxsulcaldepth = 70;           %  Maximum Sulcal Depth
    end
    if ~isfield(opts,'maxsulcallength')
        opts.maxsulcallength = 200;         % Maximum Sulcal Length
    end
    if ~isfield(opts,'ncurvp')
        opts.ncurvp = 40;     % Number of curve points
    end
    if ~isfield(opts,'verbose')
        opts.verbose = 0;     % Number of curve points
    end
end

if nargin < 3
    error('Three Inputs are needed');
    return
end
SulcSurfMat = varargin{1};
PialSurfMat = varargin{2};
HullSurfMat = varargin{3};
%% ==================== End of Input parameters  =========================%

% % % % % % load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab_single.mat');
%% =========================== Main Program ============================= %

% Computing Maximum Edge length
Surfj = Compound_Surf(SulcSurfMat);

maxdist =             sqrt(max(sum((Surfj.SurfData.vertices(Surfj.SurfData.faces(:,1),:)-Surfj.SurfData.vertices(Surfj.SurfData.faces(:,2),:)).^2,2)));
maxdist = max(maxdist,sqrt(max(sum((Surfj.SurfData.vertices(Surfj.SurfData.faces(:,1),:)-Surfj.SurfData.vertices(Surfj.SurfData.faces(:,3),:)).^2,2))));
maxdist = max(maxdist,sqrt(max(sum((Surfj.SurfData.vertices(Surfj.SurfData.faces(:,2),:)-Surfj.SurfData.vertices(Surfj.SurfData.faces(:,3),:)).^2,2))));


% Hull Surface planes
VertH = HullSurfMat.SurfData.vertices;FacesH = HullSurfMat.SurfData.faces;NormalsH = HullSurfMat.SurfData.VertexNormals;
nh = cross(VertH(FacesH(:,1),:)-VertH(FacesH(:,2),:),VertH(FacesH(:,3),:)-VertH(FacesH(:,2),:),2); DH = -1*dot(nh,VertH(FacesH(:,1),:),2);

%% ==================== Sulcus Walls Labelling ========================== %
for i = 1:length(SulcSurfMat)
    Surf2Process(i) = Sulcal_Face_Labelling(SulcSurfMat(i));
end

Surfj = Compound_Surf(Surf2Process);
Graph = surface2graph(Surfj,1);
LabNet = Label_Graph_Components(Graph);
%% ================== End of Sulcus Walls Labelling ===================== %

col = [[1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];[213 221 227;175 206 227;149 196 228; 120 186 230;87 172 231;24 146 232;6 73 119;244 207 154;244 192 117;244 179 86;244 161 43; 212 133 20;158 101 19; 113 71 12]/255];
% ----------------------------------------------------------- %
%% =============== Computing Sulcal Length ============================= %%

for i = 1:length(Surf2Process)
    
    if opts.verbose
        disp(' ');
        disp('Extracting the interception curve between the WM Hull Surface and Sulcus median mesh ...');
        tic;
    end
    [sulchullIntercep, sulchullIntNormals] = Intercept_Surface_with_Surface(HullSurfMat,  Surf2Process(i));
    if opts.verbose
        toc;
    end
    if size(sulchullIntercep,1) < 5 % At least 5 interception points are needed to construct a correct interception curve
        sulchullIntercep = '';
    end
    if ~isempty (sulchullIntercep)
        %     topCoords = fitCurveTo3DPts(sulchullIntercep(:,1:3), sulchullIntercep(maxdistPair(1),1:3), sulchullIntercep(maxdistPair(2),1:3),size(sulchullIntercep,1)*2, 0 ); % Fitting a 3D curve
        [topCoords,intNormals, intBinormals] = Compute_mean_between_walls(sulchullIntercep,sulchullIntNormals);
        
        if opts.verbose
            disp(' ');
            disp('Creating perpendicular planes to the interception curve... ');
            tic;
        end
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
        intNormals(ind,:) = [];
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
            intNormals(ind,:) = [];
        end
        
        if length(unique(topCoords(:,4))) >1
            topCoords = topCoords(:,1:3) + [rand(size(topCoords,1),1) rand(size(topCoords,1),1) rand(size(topCoords,1),1)]*eps; % Adding some noise to avoid parallel planes
            temptopCoords = topCoords;
            topCoords = fitCurveTo3DPts(temptopCoords(:,1:3), temptopCoords(1,1:3), temptopCoords(end,1:3),size(temptopCoords,1)+floor(size(temptopCoords,1)*0.2), 0 ); % Fitting a 3D curve
            
            if nargin > 1
                %            xtempNormal = interp1(temptopCoords(:,1),intNormals(:,1),topCoords(:,1));
                %            ytempNormal = interp1(temptopCoords(:,2),intNormals(:,3),topCoords(:,2));
                %            ztempNormal = interp1(temptopCoords(:,3),intNormals(:,2),topCoords(:,3));
                %            intNormals = [xtempNormal(:) ytempNormal(:) ztempNormal(:)];
                FO = fit(temptopCoords(:,1),intNormals(:,1),'smoothingspline');
                xtempNormal = FO(topCoords(:,1));
                FO = fit(temptopCoords(:,2),intNormals(:,2),'smoothingspline');
                ytempNormal = FO(topCoords(:,2));
                FO = fit(temptopCoords(:,3),intNormals(:,3),'smoothingspline');
                ztempNormal = FO(topCoords(:,3));
                
            end
        end
        %% =================== Computing Top and Bottom Lines =============== %
        [bottomlineCoords, toplineCoords] =  Compute_Sulcal_Top_Bottom_Lines(HullSurfMat,Surfj,topCoords, intNormals, opts);
        
        hold on
        plot3(bottomlineCoords(:,1),bottomlineCoords(:,2),bottomlineCoords(:,3),'-','Linewidth',5,'Color',col(i,:));
        plot3(toplineCoords(:,1),toplineCoords(:,2),toplineCoords(:,3),'-','Linewidth',5,'Color',col(i,:));
    end
end

    
    %% ============ End of Computing Top and Bottom Lines =========== %
    
        
    %% =================== Computing Depth Lines ==================== %
     
    [SulcVert] =  Compute_Sulci_Depth(Surf2Process,bottomlineCoords,toplineCoords,opts);
    Nplanes = size(bottomlineCoords,1);
    opts.ncurvp = Nplanes;
    ndepth = size(SulcVert,1)/Nplanes;
    DlinesP = [0 0];
    DlinesF = [0 0];
    for planes = 1:Nplanes
        DepthCoords = SulcVert((planes-1)*ndepth+1:planes*ndepth,:);
        
        if (planes == 1)|(planes == Nplanes)
            DlinesF = [DlinesF;[(1:(ndepth-1))' (2:ndepth)'] + ndepth*(planes-1)];
        else
            DlinesP = [DlinesP;[(1:(ndepth-1))' (2:ndepth)'] + ndepth*(planes-1)];
        end
        
        SulcMetrics.depth.Coords{planes} = DepthCoords;
        DepthProfile(:,planes) = cumsum([0 sqrt(sum(((DepthCoords(1:end-1,:) - DepthCoords(2:end,:)).^2)'))]');
        
    end
    
    for planes = 1:ndepth
        Temp = [planes:ndepth:(ndepth*(size(bottomlineCoords,1)-1)+planes)]';
        if (planes >2) | (planes < ndepth)
            DlinesP = [DlinesP;[Temp(1:end-1) Temp(2:end)]];
        end
        
        templine = SulcVert(Temp(1:end),:);
        LengthProfile(planes,:) = cumsum([0 sqrt(sum(((templine(1:end-1,:) - templine(2:end,:)).^2)'))]')';
    end
    DlinesP(1,:) = [];
    DlinesF(1,:) = [];
    
    % Sulcal Length
    
    % Computing Interception Length
    verts = topCoords(:,1:3);
    intercpLength = sqrt(sum((verts(1:end-1,:) - verts(2:end,:)).^2,2));
    
    % Computing Length Metrics
    tempLengthprof = LengthProfile(:,end);
    if opts.maxsulcallength
        ind = find(tempLengthprof < opts.maxsulcallength);
        tempLengthprof = tempLengthprof(ind);
    end
    
    
    SulcMetrics.length.measures = [sum(intercpLength); max(tempLengthprof(:)); min(tempLengthprof(:)); mode(tempLengthprof(:)) ; median(tempLengthprof(:)) ;  mean(tempLengthprof(:));  std(tempLengthprof(:))];
    SulcMetrics.length.profile = LengthProfile;
    
    
    % Computing Depth Metrics
    tempDepthprof = DepthProfile(end,:);
    if opts.maxsulcaldepth
        ind = find(tempDepthprof < opts.maxsulcaldepth);
        tempDepthprof = tempDepthprof(ind);
    end
    
    % Sulcal Depth
    SulcMetrics.depth.measures = [max(tempDepthprof); min(tempDepthprof); mode(tempDepthprof) ; median(tempDepthprof) ;  mean(tempDepthprof);  std(tempDepthprof)];
    SulcMetrics.depth.profile = DepthProfile;

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
    
    [Surfwidth, sulcalWidth] = Compute_Sulcal_Spam(PialSurfMat, Surfo, opts);
    WidthProfile = reshape(sulcalWidth,[ndepth Nplanes]);
    
    % Computing Width Metrics
    tempWidthprof = WidthProfile(:);
    if opts.maxsulcalwidth
        ind = find(tempWidthprof < opts.maxsulcalwidth);
        tempWidthprof = tempWidthprof(ind);
    end
    
    SurfL(5) = Surfwidth;
    SulcMetrics.width.measures = [max(tempWidthprof); min(tempWidthprof) ; mode(tempWidthprof) ; median(tempWidthprof)  ;  mean(tempWidthprof);  std(tempWidthprof)];
    SulcMetrics.width.profile = WidthProfile;
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
    sts = unique(sulchullIntercep(:,4));
    for i = 1:length(sts)
        ind = find(sulchullIntercep(:,4) == sts(i));
        Nintvert = length(ind);
        if i ==1
            Surft.SurfData.vertices = sulchullIntercep(ind,1:3);
            Surft.SurfData.faces = [[1:Nintvert-1]' [2:Nintvert]'];
        else
            Surft.SurfData.faces = [Surft.SurfData.faces;[[1:Nintvert-1]' [2:Nintvert]']+size(Surft.SurfData.vertices,1)];
            Surft.SurfData.vertices = [Surft.SurfData.vertices;sulchullIntercep(ind,1:3)];
        end
    end
    SurfL(6) = Surft;
    
    % Saving Mean Interception Line
    Nintvert = size(topCoords,1);
    Surft.SurfData.vertices = topCoords;
    Surft.SurfData.faces = [[1:Nintvert-1]' [2:Nintvert]'];
    SurfL(7) = Surft;
    
    %% =============== End of Creating a New Surface =============== %%
    %% ================= End of Computing Depth Lines =================== %
% else
%     % Reparametrized Surface
%     Surfo = '';
%     
%     % Curves
%     SurfL = '';
%     
%     % Sulcal Width
%     SulcMetrics.width.measures = zeros(1,6);
%     SulcMetrics.width.profile = 0;
%     % Sulcal Length
%     SulcMetrics.length.measures = zeros(1,7);
%     SulcMetrics.length.profile = 0;
%     
%     % Sulcal Depth
%     SulcMetrics.depth.measures = zeros(1,6);
%     SulcMetrics.depth.profile = 0;
% end

varargout{1} = SulcMetrics;
varargout{2} = Surfo;
varargout{3} = SurfL;
return

%% ======================= Internal Functions =========================== %



function norma = normm(M)
norma = sqrt(sum((M').^2))';
return