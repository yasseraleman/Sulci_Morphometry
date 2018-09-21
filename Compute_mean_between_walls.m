function varargout = Compute_mean_between_walls(varargin);
%
% Syntax :
%       [topCoords, interpNormals, interpBinormals] = Compute_mean_between_walls(spacCoord, spacNormals);
%
% This function computes the mean interception curve between a plane and a sulcus median mesh.
%
%
% Input Parameters:
%       spacCoord               : 3D coordinates of the interception curve.
%       spacNormals             : Normals to the interception curve.
%
% Output Parameters:
%       topCoords               : Mean interception curve.
%     interpNormals             : Normals over the mean interception curve.
%     interpBinormals           : BiNormals to the interception curve (see frenet coordinate system).
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% February 13th 2015
% Version $1.0

%% =========================== Input parameters  =========================%
% load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab_test1.mat');
spacCoord = varargin{1};
if nargin >1
    spacNormals = varargin{2};
end
%% ===================== End of Input parameters  ========================%

%% ======================== Main Program  ================================%

[maxdistPair] = Sort_3D_pointCloud(spacCoord); % Detecting extreme points for each cluster

Nc = size(maxdistPair,1); % Number of cluster
topCoords = [0 0 0 0];
interpNormals = [0 0 0];
interpBinormals = [0 0 0];
extremPairs = [0 0 0 0];
for i = 1:Nc % For each cluster
    
    % Selecting the cluster points
    indp = find(spacCoord(:,4) == i); 
    sulchullIntercepCl = spacCoord(indp,:);
    if nargin == 2
       sulchullIntNormalsCl = spacNormals(indp,:);
    end
    maxdistPairCl = maxdistPair(i,:); % Extremes
    [aTemp,bTemp] = ismember(indp,maxdistPairCl(1:2));
    exTremes(nonzeros(bTemp)) = find(aTemp);
    
    [a.b] = ismember(indp,maxdistPairCl(1:2));
    
    Npoints = size(sulchullIntercepCl,1); % Number of points in each cluster
    
    indPairs = [1:Npoints-1]; %Creating conections
    indPairs = [indPairs(:) indPairs(:)+1;indPairs(end)+1 indPairs(1)];
    
    grapH = sparse(Npoints,Npoints); % Creating Graph
    ind = sub2ind([Npoints Npoints],[indPairs(:,1) indPairs(:,2)],[indPairs(:,2) indPairs(:,1)]);
    grapH(ind) =  1;
    
    % Detecting one side curve
    [distLeft,pathLeft] = graphshortestpath(grapH,exTremes(1),exTremes(2));
    ind2del = sub2ind([Npoints Npoints],[pathLeft(1:end-1) pathLeft(2:end)],[pathLeft(2:end) pathLeft(1:end-1)]); % Deleting the conections to obtain the other side
    grapH(ind2del) =  0;

    % Detecting the other side curve
    [distRight,pathRight] = graphshortestpath(grapH,exTremes(1),exTremes(2));

    if distLeft < distRight*0.4 % Using only one of both sides
        tempPoints = sulchullIntercepCl(pathRight,1:3);
%         topCoordsCl = fitCurveTo3DPts(tempPoints, tempPoints(1,:), tempPoints(end,:),length(pathRight), 0 ); % Fitting a 3D curve
               
        topCoordsCl = interparc(0:1/(length(pathRight)-1):1,tempPoints(:,1),tempPoints(:,2),tempPoints(:,3),'linear');

        
        [T,N,B,k,t] = frenet(topCoordsCl(:,1),topCoordsCl(:,2),topCoordsCl(:,3));

        if nargin > 1
           tempNormal = sulchullIntNormalsCl(pathRight,:);
%            xtempNormal = interp1(tempPoints(:,1),tempNormal(:,1),topCoordsCl(:,1));
%            ytempNormal = interp1(tempPoints(:,2),tempNormal(:,2),topCoordsCl(:,2));
%            ztempNormal = interp1(tempPoints(:,3),tempNormal(:,3),topCoordsCl(:,3));

            FO = fit(tempPoints(:,1),tempNormal(:,1),'smoothingspline');
            xtempNormal = FO(topCoordsCl(:,1));
            FO = fit(tempPoints(:,2),tempNormal(:,2),'smoothingspline');
            ytempNormal = FO(topCoordsCl(:,2));
            FO = fit(tempPoints(:,3),tempNormal(:,3),'smoothingspline');
            ztempNormal = FO(topCoordsCl(:,3));

           interpNormalsCl = [xtempNormal(:) ytempNormal(:) ztempNormal(:)];
           tempVar = normm(interpNormalsCl);
           interpNormalsCl = interpNormalsCl./[tempVar tempVar tempVar];
        else
            interpNormalsCl = N;
        end
        intBinormalsCl = B;
    elseif distRight < distLeft*0.4 % Using only one of both sides
        tempPoints = sulchullIntercepCl(pathLeft,1:3);
%         topCoordsCl = fitCurveTo3DPts(tempPoints, tempPoints(1,:), tempPoints(end,:),length(pathLeft), 0 ); % Fitting a 3D curve
        
        topCoordsCl = interparc(0:1/(length(pathLeft)-1):1,tempPoints(:,1),tempPoints(:,2),tempPoints(:,3),'linear');

        [T,N,B,k,t] = frenet(topCoordsCl(:,1),topCoordsCl(:,2),topCoordsCl(:,3));

        if nargin > 1
            tempNormal = sulchullIntNormalsCl(pathLeft,:);
%             xtempNormal = interp1(tempPoints(:,1),tempNormal(:,1),topCoordsCl(:,1));
%             ytempNormal = interp1(tempPoints(:,2),tempNormal(:,2),topCoordsCl(:,2));
%             ztempNormal = interp1(tempPoints(:,3),tempNormal(:,3),topCoordsCl(:,3));

            FO = fit(tempPoints(:,1),tempNormal(:,1),'smoothingspline');
            xtempNormal = FO(topCoordsCl(:,1));
            FO = fit(tempPoints(:,2),tempNormal(:,2),'smoothingspline');
            ytempNormal = FO(topCoordsCl(:,2));
            FO = fit(tempPoints(:,3),tempNormal(:,3),'smoothingspline');
            ztempNormal = FO(topCoordsCl(:,3));
            
            interpNormalsCl = [xtempNormal(:) ytempNormal(:) ztempNormal(:)];
            tempVar = normm(interpNormalsCl);
            interpNormalsCl = interpNormalsCl./[tempVar tempVar tempVar];
        else
            interpNormalsCl = N;
        end
        intBinormalsCl = B;
    else % Averaging both sides
       tempPoints = sulchullIntercepCl(pathRight,1:3);
       intNpoints = max([length(pathRight) length(pathLeft)]);
%        topCoordsCl = fitCurveTo3DPts(tempPoints, tempPoints(1,:), tempPoints(end,:),intNpoints, 0 ); % Fitting a 3D curve
       
       topCoordsCl = interparc(0:1/(intNpoints-1):1,tempPoints(:,1),tempPoints(:,2),tempPoints(:,3),'linear');

       if nargin > 1
            tempNormal = sulchullIntNormalsCl(pathRight,:);
            
            FO = fit(tempPoints(:,1),tempNormal(:,1),'smoothingspline');
            xtempNormal = FO(topCoordsCl(:,1));
            FO = fit(tempPoints(:,2),tempNormal(:,2),'smoothingspline');
            ytempNormal = FO(topCoordsCl(:,2));
            FO = fit(tempPoints(:,3),tempNormal(:,3),'smoothingspline');
            ztempNormal = FO(topCoordsCl(:,3));
            
%             xtempNormal = interp1(tempPoints(:,1),tempNormal(:,1),topCoordsCl(:,1));
%             ytempNormal = interp1(tempPoints(:,2),tempNormal(:,2),topCoordsCl(:,2));
%             ztempNormal = interp1(tempPoints(:,3),tempNormal(:,3),topCoordsCl(:,3));
            interpNormalsCl = [xtempNormal(:) ytempNormal(:) ztempNormal(:)];
        end
       
       tempPoints = sulchullIntercepCl(pathLeft,1:3);
%        tempCoord = fitCurveTo3DPts(tempPoints, tempPoints(1,:), tempPoints(end,:),intNpoints, 0 );
       warning off;
       tempCoord = interparc(0:1/(intNpoints-1):1,tempPoints(:,1),tempPoints(:,2),tempPoints(:,3),'linear');
       
       
       % Binormal
       intBinormalsCl = tempCoord - topCoordsCl; 
       tempVar = normm(intBinormalsCl);
       intBinormalsCl(2:end-1,:) = intBinormalsCl(2:end-1,:)./[tempVar(2:end-1) tempVar(2:end-1) tempVar(2:end-1)];
       intBinormalsCl(1,:) = intBinormalsCl(2,:);
       intBinormalsCl(end,:) = intBinormalsCl(end-1,:);
       
       % Mean Interception curve
       topCoordsCl = topCoordsCl + tempCoord; % Fitting a 3D curve
       topCoordsCl = topCoordsCl/2;
       
       % Normals
       if nargin > 1
            tempNormal = sulchullIntNormalsCl(pathLeft,:);
            
            FO = fit(tempPoints(:,1),tempNormal(:,1),'smoothingspline');
            xtempNormal = interpNormalsCl(:,1) + FO(topCoordsCl(:,1));

            FO = fit(tempPoints(:,2),tempNormal(:,2),'smoothingspline');
            ytempNormal = interpNormalsCl(:,2) + FO(topCoordsCl(:,2));

            FO = fit(tempPoints(:,3),tempNormal(:,3),'smoothingspline');
            ztempNormal = interpNormalsCl(:,3) + FO(topCoordsCl(:,3));


            
%             xtempNormal = interpNormalsCl(:,1) + interp1(tempPoints(:,1),tempNormal(:,1),topCoordsCl(:,1));
%             ytempNormal = interpNormalsCl(:,2) + interp1(tempPoints(:,2),tempNormal(:,2),topCoordsCl(:,2));
%             ztempNormal = interpNormalsCl(:,3) + interp1(tempPoints(:,3),tempNormal(:,3),topCoordsCl(:,3));
            
            
            
            interpNormalsCl = [xtempNormal(:) ytempNormal(:) ztempNormal(:)]/2; % Normals
            tempVar = normm(interpNormalsCl);
            interpNormalsCl = interpNormalsCl./[tempVar tempVar tempVar];
        end
    end
    dCl = sum(sqrt(sum((topCoordsCl(1:end-1,:) - topCoordsCl(2:end,:)).^2,2))); % Mean Line geodesic length
    topCoords = [topCoords;topCoordsCl i*ones(size(topCoordsCl,1),1)];  
    interpBinormals = [interpBinormals;intBinormalsCl];
    if nargin == 2
        interpNormals = [interpNormals;interpNormalsCl];
    end
    indc = find(topCoords(:,4) == i); 
    extremPairs = [extremPairs;[indc(1) indc(end)]-1 i dCl];
end
topCoords(1,:)   = [];
interpBinormals(1,:) = [];
if nargin > 1
    interpNormals(1,:) = [];
end
extremPairs(1,:) = [];

if size(extremPairs,1) >1
    % Sorting lines
    tempVar = extremPairs(:,1:2)';
    tempVar = tempVar(:);
    extremCoords = topCoords(tempVar,1:3);
    
    distMat = dist(extremCoords'); % Euclidean Distance
%     [X,Y] = find(distMat(1,:) == max(distMat(1,:)));
    [X,Y] = find(distMat == max(distMat(:)));
    globExtrems = tempVar([X(1) Y(1)]);
    
    refPoint = globExtrems(1);
    refCoord = topCoords(refPoint,1:3);
    dis2refPoint = sqrt(sum((topCoords(:,1:3) - repmat(refCoord,[size(topCoords,1) 1])).^2,2));
    [~,ord] = sort(dis2refPoint);
    topCoords = topCoords(ord,:);
    if ~isempty(interpBinormals)
       interpBinormals = interpBinormals(ord,:);
    end
    if nargin > 1
        interpNormals = interpNormals(ord,:);
    end
    cont = 0;
    reordtopCoords = [0 0 0 0];
    reordinterpNormals = [0 0 0];
    reordinterpBinormals = [0 0 0];
    extremPairs = [0 0 0];
    Nprev = 0;
    while ~isempty(topCoords)
        cont = cont + 1;
        ind = find(topCoords(:,4) == topCoords(1,4));
        reordtopCoords = [reordtopCoords;[topCoords(ind,1:3) cont*ones(length(ind),1)]];
        topCoords(ind,:) = [];
        
        if ~isempty(interpBinormals)
            reordinterpBinormals = [reordinterpBinormals;interpBinormals(ind,1:3)];
            interpBinormals(ind,:) = [];
        end

        if nargin > 1
            reordinterpNormals = [reordinterpNormals;interpNormals(ind,1:3)];
            interpNormals(ind,:) = [];
        end
        extremPairs = [extremPairs; [ind(1)+Nprev ind(end)+Nprev cont]];
        Nprev = Nprev + length(ind);
    end
    topCoords = reordtopCoords(2:end,:);
    
    if ~isempty(interpBinormals)
       interpBinormals = reordinterpBinormals(2:end,:);
    end
    
    if nargin > 1
       interpNormals = reordinterpNormals(2:end,:);
    end
    extremPairs = extremPairs(2:end,:);
end
%% ===================== End of Main Program  ============================%
% Output
varargout{1} = topCoords;
if sum(interpNormals(:)) ==0
    interpNormals = '';
end
varargout{2} = interpNormals;
varargout{3} = interpBinormals;


return