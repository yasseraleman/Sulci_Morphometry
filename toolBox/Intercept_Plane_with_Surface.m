function varargout = Intercept_Plane_with_Surface(varargin);
%
% Syntax :
%      [Xline, Yline, Zline, clustId] = Intercept_Plane_with_Surface(planEq, Surf);
%
% This function cuts a surface (defined by vertices and faces as represented
% by MATLAB patches) by a plane, leading to a curve in 3D space. The
% resulting curve is represented by a set of contigous lines in the space
%
% Input Parameters:
%       planEq                  : A 4-length vector with the parameters of the plane. If plane = [A
%                                 B C D] then every 3D point P = (x,y,z) belonging to the plane satisfies
%                                 plane*[P; 1]' = A*x + B*y + C*z + D = 0
%       Surf                    : Surface structure as represented in MATLAB by patches:
%                                 Surf.SurfData.vertices
%                                 Surf.SurfData.faces
%
% Output Parameters:
%      intCurve                 : Interception curve.
%                                 
%       edgeLabels              : Interception edges labeled according to
%                                 the its endpoints labels the interception
%                                 surface.
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

% load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab_test1.mat');
if nargin <2
    error('Wrong inputs');
    return;
end
planEq = varargin{1};
Surf = varargin{2};
fv = Surf.SurfData;

warning off %#ok
Xline = []; Yline = []; Zline = [];
oo = ones(size(fv.vertices,1),1);
maxdist = sqrt(max(sum((fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,2),:)).^2,2)));
maxdist = max(maxdist,sqrt(max(sum((fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,3),:)).^2,2))));
maxdist = max(maxdist,sqrt(max(sum((fv.vertices(fv.faces(:,2),:)-fv.vertices(fv.faces(:,3),:)).^2,2))));
vertx = find(abs(dot([fv.vertices oo],planEq(oo,:),2)/norm(planEq(1:3)))<maxdist);
indf = ismember(fv.faces,vertx);
[rindf,cindf] = find(indf); %#ok
rindf = unique(rindf);
Nf = length(rindf);
temp = [1 2;2 3; 3 1];
faceAndEdges = [0 0 0];
for i = 1:Nf
    verts = fv.vertices(fv.faces(rindf(i),:),:);
    intPoints = Intercept_Plane_with_LinesSegment(planEq, verts(temp(:,1),1:3),verts(temp(:,2),1:3));
    if ~isempty(intPoints)&(size(intPoints,1)==2)
         Xline = [Xline intPoints(:,1)]; %#ok
         Yline = [Yline intPoints(:,2)]; %#ok
         Zline = [Zline intPoints(:,3)]; %#ok
         faceAndEdges  = [faceAndEdges; [[fv.faces(rindf(i),temp(intPoints(1,4),:)); fv.faces(rindf(i),temp(intPoints(2,4),:))]  [rindf(i);rindf(i)]]];
    end
end
faceAndEdges(1,:) = [];

temp = sort(faceAndEdges(:,1:2)')';
tempXline = Xline(:);
tempYline = Yline(:);
tempZline = Zline(:);
[finintEdges,tempindex,indexes] = unique(temp(:,1:2),'rows'); % Final Interception edges
tempXline = tempXline(tempindex);
tempYline = tempYline(tempindex);
tempZline = tempZline(tempindex);

origEdges = finintEdges(:,1:2);
intFaces = fv.faces(rindf,:);

allEdges = [intFaces(:,1) intFaces(:,2); intFaces(:,2) intFaces(:,3);intFaces(:,1) intFaces(:,3)] ; % Edges from the intersection faces

allEdges = sort(allEdges')';  % Unifying edges to remove repeated edges
allEdges = [allEdges repmat(rindf,[3 1])]; % Labelling each edge according to its face.

% Detecting the faces that connects each pair of edges
Nv = size(origEdges,1);
sharedFaces = zeros(Nv,2);


for i = 1:Nv;
    ind = find(sum(ismember(intFaces,origEdges(i,:)),2)==2);
    sharedFaces(i,:) = ind(:)';
end


definitEdges = full(Reorganize_Intercep_Between_Surfaces(Surf,finintEdges));
[~,neworder] = ismember(definitEdges(:,1:2),finintEdges(:,1:2),'rows');
tempXline = tempXline(neworder);
tempYline = tempYline(neworder);
tempZline = tempZline(neworder);
edgeId = definitEdges(:,3);

intCurve = [tempXline tempYline tempZline edgeId];
% % 
% % idClusters = nonzeros(unique(edgeId(:))); % Number of possible
% % Nc = length(idClusters); 
% % for i  = 1:Nc
% %     ind = find(edgeId == idClusters(i)); % Selecting the cluster
% %     if i == 1
% %         Xline = [tempXline(ind(1:end-1)) tempXline(ind(2:end))]';
% %         Yline = [tempYline(ind(1:end-1)) tempYline(ind(2:end))]';
% %         Zline = [tempZline(ind(1:end-1)) tempZline(ind(2:end))]';
% %         clustId = ones(1,length(ind)-1);
% %     else
% %         Xline = [Xline [tempXline(ind(1:end-1)) tempXline(ind(2:end))]'];
% %         Yline = [Yline [tempYline(ind(1:end-1)) tempYline(ind(2:end))]'];
% %         Zline = [Zline [tempZline(ind(1:end-1)) tempZline(ind(2:end))]'];
% %         clustId = [clustId ones(1,length(ind)-1)*i];
% %     end
% % end
if isfield(Surf,'Is')
    edgeLabels = Surf.Is(definitEdges(:,1:2));
else
    edgeLabels = '';
end
%% ================= End of Reordering the interception edges =========== %
%% =========================== End of Main Program ====================== %
varargout{1} = intCurve;
varargout{2} = edgeLabels;


warning on %#ok
return