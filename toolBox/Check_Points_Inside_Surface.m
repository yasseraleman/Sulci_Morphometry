function varargout = Check_Points_Inside_Surface(varargin)
%
% Syntax :
% inVec = Check_Points_Inside_Surface(Surf,pCoords);
%
% This function tests if points are inside a 3D triangulated (faces/vertices) surface
%
% Input Parameters:
%   Surf        : Surfaces files.
%   pCoords     : Points coordinates
%
% Output Parameters:
%   inVec       : Boolean vector containing 1s for points inside the
%                 surface and 0 otherwise.
%
% Related references:
%
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 30th 2006
% Version $1.0

%% ====================== Checking input parameters ===================== %
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    Surf = varargin{1};
    Surf = Surface_Checking(Surf);
    Surf = Surf(1);
end
if nargin <2
    error('You must enter points coordinates');
end
pCoords = varargin{2};
%% =============== End of Checking input parameters ===================== %

%% ======================== Main Program  =============================== %
% Surface vertices and faces
VertP  = Surf.SurfData.vertices;
FacesP = Surf.SurfData.faces;
np = cross(VertP(FacesP(:,1),:)-VertP(FacesP(:,2),:),VertP(FacesP(:,3),:)-VertP(FacesP(:,2),:),2); Dp = -1*dot(np,VertP(FacesP(:,1),:),2);

% Random Direction
% rDirec = rand(1,3);rDirec = rDirec./sqrt(rDirec*rDirec');
rDirec = [0 0 1];
rDirec = repmat(rDirec,[size(pCoords,1) 1]);

% Creating line coordinates
P1 = pCoords; P2 = pCoords-rDirec;
Npoints = size(P1,1);

inVec = zeros(size(pCoords,1),1);
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
    
    if ~isempty(inte)
        orientsign = sign((dot(inte -repmat(pCoords(po,:),[size(inte,1) 1]),repmat(rDirec(po,:),[size(inte,1) 1]),2))); % Detecting Orientation
        ind = find(orientsign >0);
        if mod(length(ind),2)
            inVec(po,1) = 1;
        end
    end
end
%% ================== End of Main Program  ============================== %
% Outputs
varargout{1} = inVec;
return;