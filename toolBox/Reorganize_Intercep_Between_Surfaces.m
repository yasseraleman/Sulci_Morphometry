function varargout = Reorganize_Intercep_Between_Surfaces(varargin)
%
% Syntax :
%     [reordintEdges] = Reorganize_Intercep_Between_Surfaces(intSurf, intEdges);
%
% This script reorganizes the interception points between surfaces
%
% Input Parameters:
%       intSurf                 : Sulcus Surface in matlab format
%       intEdges                : Interception Edges
%
% Output Parameters:
%       reordintEdges          : Reordered interception edges. The third
%                                row shows the cluster Id for each
%                                interception edge
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin < 2
    error('There are missing input variables');
    return;
end

%% =========================== Main Program ============================= %

intSurf  = varargin{1}; % Interception Surface
intEdges = varargin{2}; % Interception Edges
origEdges = intEdges;
intFaces = intSurf.SurfData.faces;

allEdges = [intFaces(:,1) intFaces(:,2); intFaces(:,2) intFaces(:,3);intFaces(:,1) intFaces(:,3)] ; % Edges from the intersection faces

allEdges = sort(allEdges')';  % Unifying edges to remove repeated edges
allEdges = [allEdges repmat([1:size(intSurf.SurfData.faces,1)]',[3 1])]; % Labelling each edge according to its face.


% Detecting the faces that connects each pair of edges
Nv = size(intEdges,1);
% Creating Graph. Edges that belongs to the same face are connected
Graph = sparse(Nv,Nv);
cont = 0;
for i = 1:Nv;
    ind = find(sum(ismember(intSurf.SurfData.faces,intEdges(i,:)),2)==2);
    for j = 1:length(ind)
        tempface = intSurf.SurfData.faces(ind(j),:);
        tempEdges = [tempface(:,1) tempface(:,2); tempface(:,2) tempface(:,3);tempface(:,1) tempface(:,3)];
        tempEdges = sort(tempEdges')';
        X = find(ismember(intEdges,tempEdges,'rows'));
        Graph(X(1),X(2)) = 1;
        Graph(X(2),X(1)) = 1;
    end
end


% Detecting if the sulcus is cutted in two pieces by the hull surface
LabNet = Label_Graph_Components(Graph);

idClusters = nonzeros(unique(LabNet(:))); % Number of possible 
Nc = length(idClusters); % Number of clusters
edgeId = max(LabNet);% Labelling each edge according to its cluster number
edgeId = edgeId(:);


%% ================= Reordering the interception edges ================== %
neworder = 0;
cont = 0;
for i  = 1:Nc
    ind = find(edgeId == idClusters(i)); % Selecting the cluster
    subGraph = Graph(ind,ind); % Creating a subgraph for each cluster
    
    if graphisdag(subGraph)
        start_point = Xrem(1); % Selecting the start and the endpoint
        end_point = Yrem(1);
        subGraph(start_point,end_point) = 0; % Disconecting the cluster
        subGraph(end_point,start_point) = 0;
        [DIST, PATHS]=graphshortestpath(subGraph,start_point,end_point); % Detecting the order between points
        PATHS = [PATHS(end) PATHS(1:end-1)];
        PATHS = PATHS(:);
        neworder = [neworder;ind(PATHS(:))]; % New order
    else
        cont = cont + 1;
        distances=graphallshortestpaths(subGraph);
        distances(find(distances==Inf))=0;
        [start_point,end_point] = find(distances == max(distances(:)));
        [DIST, PATHS]=graphshortestpath(subGraph,start_point(1),end_point(1)); % Detecting the order between points
        
        if cont == 1
            extPairs(cont,:) = [1 length(PATHS) idClusters(i)];
        else
            extPairs(cont,:) = [length(neworder) length(neworder)+length(PATHS)-1 idClusters(i)];
        end
        neworder = [neworder;ind(PATHS(:))]; % New order
    end
end
neworder(1) = [];
reordintEdges = [intEdges(neworder,:) edgeId(neworder,:)];


if exist('extPairs','var')
    if size(extPairs,1) >1
        
        Graph = surface2graph(intSurf);
        distances=graphallshortestpaths(Graph);
        distances(find(distances==Inf))=0;
        [extremes,fillZero] = unique(extPairs(:,1:2));
        fillZero = reshape(fillZero,[2 size(extPairs,1)])';
        Nlen = length(extremes);
        indfillZero = sub2ind([Nlen Nlen],[fillZero(:,1);fillZero(:,2)],[fillZero(:,2);fillZero(:,1)]);
        temDist = distances(reordintEdges(extremes),reordintEdges(extremes));
        [Xextr,Yextr] = find(temDist == max(temDist(:)));
        
        
        [X,Y] = meshgrid(extremes,extremes);X = X(:);Y = Y(:);
        tempOrd = unique(sort([X Y]')','rows');
        ord2rem = find(ismember(tempOrd,[[extremes extremes];extPairs(:,1:2)],'rows'));
        tempOrd(ord2rem,:) =[];
        [a,b] = ismember(tempOrd,extremes);
        ind = sub2ind(size(temDist),b(:,1),b(:,2));
        
        distVect = [extremes(b) temDist(ind)];
        
        tempexT = extPairs(:,1:2);
        tempVal = ismember(tempexT,extremes(Xextr(1)));
        neighPos = find(tempVal - repmat(sum(tempVal,2),[1 2]) == -1);
        stNeigh = tempexT(neighPos);
        reordExtremes = [extremes(Xextr(1)) stNeigh 1];
        tempexT = [tempexT;flipdim(tempexT,2)];
        for i = 1:Nc-1
            oldstNeigh = stNeigh;
            [X,Y] = find(distVect(:,1:2) == oldstNeigh);
            [~,pos] = min(distVect(X,3));
            stNeight = distVect(X(pos),1:2);
            stNeight(stNeight == oldstNeigh) = [];
            
            [X,Y] = find(tempexT(:,1) == stNeight);
            stNeigh = tempexT(X,2);
            reordExtremes = [reordExtremes;stNeight stNeigh i+1];
            %         extPairs == Xextr
        end
        %
        flipBool = (reordExtremes(:,1) - reordExtremes(:,2));
        definitEdges(1,:) = [0 0 0];
        for i = 1:Nc
            partExtremes = reordExtremes(i,:);
            indrow = find(sum(ismember(extPairs(:,1:2),partExtremes(:)),2)== 2);
            clustId = extPairs(indrow,3);
            IntercEdgesTemp = reordintEdges(find(reordintEdges(:,3) == clustId),:);
            if flipBool(i)>0
                IntercEdgesTemp = flipdim(IntercEdgesTemp,1);
            end
            definitEdges = [definitEdges;IntercEdgesTemp(:,1:2) ones(size(IntercEdgesTemp,1),1)*i];
        end
        definitEdges(1,:) = [];
        reordintEdges= definitEdges;
    end
end
%% =========================== End of Main Program ====================== %

varargout{1} = reordintEdges;
return;