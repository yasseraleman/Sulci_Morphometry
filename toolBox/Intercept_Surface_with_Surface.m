function varargout = Intercept_Surface_with_Surface(varargin);
%
% Syntax :
%       [intCurve, edgeLabels] = Intercept_Surface_with_Surface(estSurface,  intSurface);
%
% This function computes the interception line between two surfaces
%
%
% Input Parameters:
%       estSurface              : Static Surface
%       intSurface              : Intercepting Surface.
%
% Output Parameters:
%      intCurve                 : Interception Curve.
%     edgeLabels                : Interception edges labeled according to
%                                 the its endpoints labels the interception
%                                 surface.
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% February 13th 2015
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin < 2
    error('Three Inputs are needed');
    return
end
estSurface = varargin{1};
intSurface = varargin{2};

% load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab_test2.mat');
%% =========================== Main Program ============================= %
% Surface planes
if ~isfield(estSurface.SurfData,'VertexNormals')
    estSurface = Compute_Surface_Normals(estSurface);
end
if ~isfield(intSurface.SurfData,'VertexNormals')
    intSurface = Compute_Surface_Normals(estSurface);
end


VertH = estSurface.SurfData.vertices;FacesH = estSurface.SurfData.faces;NormalsH = estSurface.SurfData.VertexNormals;
nh = cross(VertH(FacesH(:,1),:)-VertH(FacesH(:,2),:),VertH(FacesH(:,3),:)-VertH(FacesH(:,2),:),2); DH = -1*dot(nh,VertH(FacesH(:,1),:),2);
% in = Check_Points_Inside_Surface(estSurface,intSurface.SurfData.vertices);
in = inpolyhedron(estSurface.SurfData,intSurface.SurfData.vertices);
indin = find(in == 1); % Vertices inside hull surface
indout = find(in == 0); % Vertices outside hull surface

intlogFaces = in(intSurface.SurfData.faces); % Detecting Interception Faces from the sulcus mesh

temp = sum(intlogFaces,2); % The faces with points in both sides of the hull surface will have a sum value minor than 3 and different from 0
indintFaces = find((temp < 3)&(temp ~=0)); % Interception Faces indexes

if ~isempty(indintFaces)
    intFaces = intSurface.SurfData.faces(indintFaces,:); % Interception Faces
    
    intfacesEdges = [intFaces(:,1) intFaces(:,2); intFaces(:,2) intFaces(:,3);intFaces(:,1) intFaces(:,3)] ; % Edges from the intersection faces
    
    temp = sort(intfacesEdges')';  % Unifying edges to remove repeated edges
    [finintEdges,tempindex,indexes] = unique(temp(:,1:2),'rows'); % Final Interception edges
    ind2keep = find(sum(in(finintEdges),2)==1);
    finintEdges = finintEdges(ind2keep,:);
%     definitEdges = Reorganize_Edges(intSurface,finintEdges);

    definitEdges = full(Reorganize_Intercep_Between_Surfaces(intSurface,finintEdges));
        
    P1 = intSurface.SurfData.vertices(definitEdges(:,1),:);
    P2 = intSurface.SurfData.vertices(definitEdges(:,2),:); % Creating Lines from points 1 and 2 (Edge 12 of the triangle)
    minvals = sqrt(sum((P1-P2).^2,2));
    
    % Number of edges
    Npoints = length(P1);
  
    intCurve = [0 0 0 0 0];
    for po =1:Npoints
        Num = dot(nh,repmat(P1(po,:),[size(nh,1) 1]),2)+DH; Den = dot(nh,repmat(P2(po,:)-P1(po,:),[size(nh,1) 1]),2); % Creating Lines
        
        t = -1*Num./Den;clear Num Den; % Lines parameters
        xint = single(repmat(P1(po,1)',[size(FacesH,1) 1])+t.*(repmat((P2(po,1)-P1(po,1))',[size(FacesH,1) 1]))); % Line parametrization
        yint = single(repmat(P1(po,2)',[size(FacesH,1) 1])+t.*(repmat((P2(po,2)-P1(po,2))',[size(FacesH,1) 1])));
        zint = single(repmat(P1(po,3)',[size(FacesH,1) 1])+t.*(repmat((P2(po,3)-P1(po,3))',[size(FacesH,1) 1])));clear t;
        
        PpP1 =  VertH(FacesH(:,1),:)-[xint yint zint];
        PpP2 =  VertH(FacesH(:,2),:)-[xint yint zint];
        PpP3 =  VertH(FacesH(:,3),:)-[xint yint zint];
        angP2p3 = (acos(dot(PpP2,PpP3,2)./(normm(PpP2).*normm(PpP3))))*180/pi; % Angles between each face point and the interception point inpial surface
        angP1p3 = (acos(dot(PpP1,PpP3,2)./(normm(PpP1).*normm(PpP3))))*180/pi;
        angP1p2 = (acos(dot(PpP1,PpP2,2)./(normm(PpP1).*normm(PpP2))))*180/pi;
        
        ind = find(round(angP2p3+angP1p3+angP1p2) == 360);
        [inte,indord] = unique([xint(ind) yint(ind) zint(ind)],'rows');ind = ind(indord); % Line interception with hull surface. ind is the intercepted face over the hull surface
        if ~isempty(inte)
            dx = (inte(:,1) - P1(po,1)).^2;
            dy = (inte(:,2) - P1(po,2)).^2;
            dz = (inte(:,3) - P1(po,3)).^2;
            dt = sqrt(dx + dy + dz);
            
            [~,pos] = min(dt);
            if dt(pos) <=minvals(po)
                intCurve = [intCurve;inte(pos,:) definitEdges(po,3) ind(pos)];
%                 plot3(P2(po,1),P2(po,2),P2(po,3),'.g','Linewidth',2,'Markersize',20);
%                 plot3(P1(po,1),P1(po,2),P1(po,3),'.r','Linewidth',2,'Markersize',20);
%                 line([P1(po,1) P2(po,1)]',[P1(po,2) P2(po,2)]',[P1(po,3) P2(po,3)]','Color',[0 0 1])
%                 plot3(inte(pos,1),inte(pos,2),inte(pos,3),'.w','Linewidth',2,'Markersize',40);
            end
        end
    end
    intCurve(1,:) = [];
    if isfield(intSurface,'Is')
        edgeLabels = intSurface.Is(definitEdges(:,1:2));
    else
        edgeLabels = '';
    end
else
    disp('These Surface do not intercept');
    intCurve = ''; 
    edgeLabels = '';
end
varargout{1} = intCurve;
varargout{2} = edgeLabels;
return