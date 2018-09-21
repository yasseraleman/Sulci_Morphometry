function varargout = Sulcal_Face_Labelling(varargin);
%
% Syntax :
% Surfj = Sulcal_Correction(SulcSurfMat,opts);
%
% This script removes anastomotic sulci from a sulci mesh. It also find
% posible sulci clusters(according to euclidean distance) and merge surfaces from the same cluster in a
% single surface. Besides, it estimate sulci surface boundary according to
% its curvature.
%
% Input Parameters:
%       SulcSurfMat             : Sulci Surface in matlab format
%       opts                    : Options
%
% Output Parameters:
%      Surfj                    : Corrected surface
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin == 0
    error('Please enter a correct Surface struct');
    return;
end

SulcSurfMat = varargin{1};
if nargin == 1
    opts.angthreshold = 40; % Angle Threshold
    opts.perc = 50;
end
if nargin == 2
    opts = varargin{2};
end
%% ====================== Identifying sulci boundaries ===================%
% disp(' ');
% disp('Identifying sulci boundaries ... ');
% tic;



for k = 1:length(SulcSurfMat)
    
     %  ------- Creating Neighoborhoods
    if ~isfield(SulcSurfMat(k),'Tri')
        Npoints = size(SulcSurfMat(k).SurfData.vertices,1);
        Nfaces = size(SulcSurfMat(k).SurfData.faces,1);
        [Tri] = Vert_Neib(double(SulcSurfMat(k).SurfData.faces),Npoints,Nfaces);
        Temp = sum(Tri);
        Tri(:,Temp==0) = [];
        SulcSurfMat(k).Tri = Tri;
    end
    
    %  ------- Smoothing the patch to obtaing strong curvatures in sulcus borders
    FV2=smoothpatch(SulcSurfMat(k).SurfData,1,50);
    
    %  ------- Computing Sulcus curvature
    [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(FV2.vertices,FV2.faces);
    temp = Compute_Surface_Normals(SulcSurfMat(k));
    
    %[tmp,normal] = compute_normal(Surfr.SurfData.vertices,Surfr.SurfData.faces);
    SulcSurfMat(k).SurfData.VertexNormals = temp.SurfData.VertexNormals;
    Normal = SulcSurfMat(k).SurfData.VertexNormals;
    
    %  ------- Separating sulcus border from the sulcus walls
    [Trip] = Vert_Neibp(double(SulcSurfMat(k).SurfData.faces),size(SulcSurfMat(k).SurfData.vertices,1),size(SulcSurfMat(k).SurfData.faces,1));
    Temp = sum(Trip);
    Trip(:,Temp==0) = [];
    temp = Trip(:,3:end);
    Nnonzeros = sum(logical(temp),2);
    indz = find(temp == 0);
    temp(indz) = 1;
    Coord = [repmat(Trip(:,1),[size(temp,2) 1]) temp(:)]; % Creating Pairs of Neighbors
    a = acos(dot(Normal(Coord(:,1),:),Normal(Coord(:,2),:),2))*180/pi; % Angle between the Normals of each pair of Neighbors
    
    b = reshape(a,[size(Trip,1) size(Trip,2)-2]);
    b(indz) = 0;
    b = (sum(b,2))./Nnonzeros; % Creating the mean of all neighborhood
    c = kmeans(b,2);
    temp = accumarray(nonzeros(c),nonzeros(c)*0+1);
    
    ind1 = find(c == 1);
    ind2 = find(c == 2);
    if b(ind1(1)) < b(ind2(1))
        pos = 1;
        SulcSurfMat(k).Is = zeros(size(SulcSurfMat(k).SurfData.vertices,1),1);
        loc = find(c ~= pos );
        SulcSurfMat(k).Is(loc) = 1;
        set(0,'RecursionLimit',1000);
        [labid] = Recur_Corr(SulcSurfMat(k),0,zeros(size(SulcSurfMat(k).Is)),1);
        if max(labid) > 2
            temp = accumarray(nonzeros(labid),nonzeros(labid)*0+1);
            [v,ord] = sort(temp,'descend');
            ind2zero = find(ismember(labid,ord(3:end)));
            c(ind2zero) = 0;
            ind = find(c);
            SulcSurfMat(k).Is = c;
            [SulcSurfMat(k)] = Surf_Corr(SulcSurfMat(k));
            SulcSurfMat(k).Is(ind) = c(ind);
            c = SulcSurfMat(k).Is;
        end
    elseif b(ind1(1)) > b(ind2(1))
        pos = 2;
        SulcSurfMat(k).Is = zeros(size(SulcSurfMat(k).SurfData.vertices,1),1);
        loc = find(c ~= pos );
        SulcSurfMat(k).Is(loc) = 1;
        set(0,'RecursionLimit',1000);
        [labid] = Recur_Corr(SulcSurfMat(k),0,zeros(size(SulcSurfMat(k).Is)),1);
        if max(labid) > 2
            temp = accumarray(nonzeros(labid),nonzeros(labid)*0+1);
            [v,ord] = sort(temp,'descend');
            ind2zero = find(ismember(labid,ord(3:end)));
            c(ind2zero) = 0;
            ind = find(c);
            SulcSurfMat(k).Is = c;
            [SulcSurfMat(k)] = Surf_Corr(SulcSurfMat(k));
            SulcSurfMat(k).Is(ind) = c(ind);
            c = SulcSurfMat(k).Is;
        end
    end
    
    ind1 = find(c == 1);
    ind2 = find(c == 2);
    if b(ind1(1)) < b(ind2(1))
        pos = 1;
    elseif b(ind1(1)) > b(ind2(1))
        pos = 2;
    end
    loc = find(c ~= pos );
    SulcSurfMat(k).Is = zeros(size(SulcSurfMat(k).SurfData.vertices,1),1);
    SulcSurfMat(k).Is(loc) = 1;
    %% ========================= Labelling Sulci Faces ===================%
    
    set(0,'RecursionLimit',1000);
    [labid] = Recur_Corr(SulcSurfMat(k),0,zeros(size(SulcSurfMat(k).Is)),1);
    SulcSurfMat(k).Is = labid;
    indzeros = find(SulcSurfMat(k).Is == 0);
    SulcSurfMat(k).Is(indzeros) = max(SulcSurfMat(k).Is) + 1;
    [SulcSurfMat(k)] = Surf_Corr(SulcSurfMat(k));
    indzeros = find(SulcSurfMat(k).Is == max(SulcSurfMat(k).Is));
    SulcSurfMat(k).Is(indzeros) = 0;

    
    % Watershed
    
    temp = Trip(:,3:end);
    indz = find(temp == 0);
    temp(indz) = 1;
    labid = SulcSurfMat(k).Is;
    Labels_Mat = labid(temp);
    Labels_Mat(indz) = 0;
    A = Indexing_Neights_V2(Trip);
    T = A(:,2:end);
    ind = find(T);
    T(ind) = T(ind) - 2*size(Trip,1);
    A(:,2:end) = T;
    indn = find(labid == 0);
    curvmap = Cmean;
    while ~isempty(indn)
        Neigh = find(sum(Labels_Mat,2) ~= 0);
        Neigh(ismember(Neigh,find(labid~=0))) = [];
        [distmin,indgrow] = min(curvmap(Neigh));
        vert_to_grow = Neigh(indgrow);
        labid(vert_to_grow) = mode(nonzeros(Labels_Mat(vert_to_grow,:)));
        L = A(vert_to_grow,1);
        Labels_Mat(A(vert_to_grow,2:L+1)) =  labid(vert_to_grow);
        indn = find(labid == 0);
    end
    
    
    if max(labid)>2
        ind2rem = find(labid > 2);
        labid(ind2rem) = 0;
        temp = Trip(:,3:end);
        indz = find(temp == 0);
        temp(indz) = 1;
        Labels_Mat = labid(temp);
        Labels_Mat(indz) = 0;
        A = Indexing_Neights_V2(Trip);
        T = A(:,2:end);
        ind = find(T);
        T(ind) = T(ind) - 2*size(Trip,1);
        A(:,2:end) = T;
        indn = find(labid == 0);
        curvmap = Cmean;
        while ~isempty(indn)
            Neigh = find(sum(Labels_Mat,2) ~= 0);
            Neigh(ismember(Neigh,find(labid~=0))) = [];
            [distmin,indgrow] = min(curvmap(Neigh));
            vert_to_grow = Neigh(indgrow);
            labid(vert_to_grow) = mode(nonzeros(Labels_Mat(vert_to_grow,:)));
            L = A(vert_to_grow,1);
            Labels_Mat(A(vert_to_grow,2:L+1)) =  labid(vert_to_grow);
            indn = find(labid == 0);
        end
        
        
    end
    SulcSurfMat(k).Is = labid;
    if max(labid) == 1
        SulcSurfMat(k) = Sulcal_Face_Labelling_kmeans3(SulcSurfMat(k));
    end
    
    %% ===================== End of Labelling Sulci Walls ================%
end
%toc;
varargout{1} = SulcSurfMat;
%% ====================== End of Identifying sulci boundaries ============%
return;


function varargout = Sulcal_Face_Labelling_kmeans3(varargin);
%
% Syntax :
% Surfj = Sulcal_Correction(SulcSurfMat,opts);
%
% This script removes anastomotic sulci from a sulci mesh. It also find
% posible sulci clusters(according to euclidean distance) and merge surfaces from the same cluster in a
% single surface. Besides, it estimate sulci surface boundary according to
% its curvature.
%
% Input Parameters:
%       SulcSurfMat             : Sulci Surface in matlab format
%       opts                    : Options
%
% Output Parameters:
%      Surfj                    : Corrected surface
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin == 0
    error('Please enter a correct Surface struct');
    return;
end

SulcSurfMat = varargin{1};
if nargin == 1
    opts.angthreshold = 40; % Angle Threshold
    opts.perc = 50;
end
if nargin == 2
    opts = varargin{2};
end


%  ------- Smoothing the patch to obtaing strong curvatures in sulcus borders
FV2=smoothpatch(SulcSurfMat.SurfData,1,50);

%  ------- Computing Sulcus curvature
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(FV2.vertices,FV2.faces);
f = figure('Visible','off');
temp = patch(FV2);

%[tmp,normal] = compute_normal(Surfr.SurfData.vertices,Surfr.SurfData.faces);
SulcSurfMat.SurfData.VertexNormals = get(temp,'VertexNormals');
norma = normm(SulcSurfMat.SurfData.VertexNormals);
SulcSurfMat.SurfData.VertexNormals = SulcSurfMat.SurfData.VertexNormals./([norma norma norma ]+eps);
Normal = SulcSurfMat.SurfData.VertexNormals;
close(f);

%  ------- Separating sulcus border from the sulcus walls
[Trip] = Vert_Neibp(double(SulcSurfMat.SurfData.faces),size(SulcSurfMat.SurfData.vertices,1),size(SulcSurfMat.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
Nnonzeros = sum(logical(temp),2);
indz = find(temp == 0);
temp(indz) = 1;
Coord = [repmat(Trip(:,1),[size(temp,2) 1]) temp(:)]; % Creating Pairs of Neighbors
a = acos(dot(Normal(Coord(:,1),:),Normal(Coord(:,2),:),2))*180/pi; % Angle between the Normals of each pair of Neighbors

b = reshape(a,[size(Trip,1) size(Trip,2)-2]);
b(indz) = 0;
b = (sum(b,2))./Nnonzeros; % Creating the mean of all neighborhood
c = kmeans(b,3);
temp = accumarray(nonzeros(c),nonzeros(c)*0+1);
[minval,ind0]  = min(temp);
[postemp]  = find(temp > minval);


ind1 = find(c == postemp(1));
ind2 = find(c == postemp(2));
if b(ind1(1)) < b(ind2(1))
    pos = 1;
    SulcSurfMat.Is = zeros(size(SulcSurfMat.SurfData.vertices,1),1);
    loc = find(c ~= pos );
    SulcSurfMat.Is(loc) = 1;
    set(0,'RecursionLimit',1000);
    [labid] = Recur_Corr(SulcSurfMat,0,zeros(size(SulcSurfMat.Is)),1);
    if max(labid) > 2
        temp = accumarray(nonzeros(labid),nonzeros(labid)*0+1);
        [v,ord] = sort(temp,'descend');
        ind2zero = find(ismember(labid,ord(3:end)));
        c(ind2zero) = 0;
        ind = find(c);
        SulcSurfMat.Is = c;
        [SulcSurfMat] = Surf_Corr(SulcSurfMat);
        SulcSurfMat.Is(ind) = c(ind);
        c = SulcSurfMat.Is;
    end
elseif b(ind1(1)) > b(ind2(1))
    pos = 2;
    SulcSurfMat.Is = zeros(size(SulcSurfMat.SurfData.vertices,1),1);
    loc = find(c ~= pos );
    SulcSurfMat.Is(loc) = 1;
    set(0,'RecursionLimit',1000);
    [labid] = Recur_Corr(SulcSurfMat,0,zeros(size(SulcSurfMat.Is)),1);
    if max(labid) > 2
        temp = accumarray(nonzeros(labid),nonzeros(labid)*0+1);
        [v,ord] = sort(temp,'descend');
        ind2zero = find(ismember(labid,ord(3:end)));
        c(ind2zero) = 0;
        ind = find(c);
        SulcSurfMat.Is = c;
        [SulcSurfMat] = Surf_Corr(SulcSurfMat);
        SulcSurfMat.Is(ind) = c(ind);
        c = SulcSurfMat.Is;
    end
end

temp = accumarray(nonzeros(c),nonzeros(c)*0+1);
[minval,ind0]  = min(temp);
[postemp]  = find(temp > minval);

ind1 = find(c == postemp(1));
ind2 = find(c == postemp(2));
if b(ind1(1)) < b(ind2(1))
    pos = postemp(1);
elseif b(ind1(1)) > b(ind2(1))
    pos =  postemp(2);
end
loc = find(c ~= pos );
SulcSurfMat.Is = zeros(size(SulcSurfMat.SurfData.vertices,1),1);
SulcSurfMat.Is(loc) = 1;
%% ========================= Labelling Sulci Faces ===================%

set(0,'RecursionLimit',1000);
[labid] = Recur_Corr(SulcSurfMat,0,zeros(size(SulcSurfMat.Is)),1);
SulcSurfMat.Is = labid;
indzeros = find(SulcSurfMat.Is == 0);
SulcSurfMat.Is(indzeros) = max(SulcSurfMat.Is) + 1;
[SulcSurfMat] = Surf_Corr(SulcSurfMat);
indzeros = find(SulcSurfMat.Is == max(SulcSurfMat.Is));
SulcSurfMat.Is(indzeros) = 0;


% Watershed

temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;
labid = SulcSurfMat.Is;
Labels_Mat = labid(temp);
Labels_Mat(indz) = 0;
A = Indexing_Neights_V2(Trip);
T = A(:,2:end);
ind = find(T);
T(ind) = T(ind) - 2*size(Trip,1);
A(:,2:end) = T;
indn = find(labid == 0);
curvmap = Cmean;
while ~isempty(indn)
    Neigh = find(sum(Labels_Mat,2) ~= 0);
    Neigh(ismember(Neigh,find(labid~=0))) = [];
    [distmin,indgrow] = min(curvmap(Neigh));
    vert_to_grow = Neigh(indgrow);
    labid(vert_to_grow) = mode(nonzeros(Labels_Mat(vert_to_grow,:)));
    L = A(vert_to_grow,1);
    Labels_Mat(A(vert_to_grow,2:L+1)) =  labid(vert_to_grow);
    indn = find(labid == 0);
end


if max(labid)>2
    ind2rem = find(labid > 2);
    labid(ind2rem) = 0;
    temp = Trip(:,3:end);
    indz = find(temp == 0);
    temp(indz) = 1;
    Labels_Mat = labid(temp);
    Labels_Mat(indz) = 0;
    A = Indexing_Neights_V2(Trip);
    T = A(:,2:end);
    ind = find(T);
    T(ind) = T(ind) - 2*size(Trip,1);
    A(:,2:end) = T;
    indn = find(labid == 0);
    curvmap = Cmean;
    while ~isempty(indn)
        Neigh = find(sum(Labels_Mat,2) ~= 0);
        Neigh(ismember(Neigh,find(labid~=0))) = [];
        [distmin,indgrow] = min(curvmap(Neigh));
        vert_to_grow = Neigh(indgrow);
        labid(vert_to_grow) = mode(nonzeros(Labels_Mat(vert_to_grow,:)));
        L = A(vert_to_grow,1);
        Labels_Mat(A(vert_to_grow,2:L+1)) =  labid(vert_to_grow);
        indn = find(labid == 0);
    end
    
    
end
SulcSurfMat.Is = labid;


%% ===================== End of Labelling Sulci Walls ================%
varargout{1} = SulcSurfMat;
%% ====================== End of Identifying sulci boundaries ============%
return;




