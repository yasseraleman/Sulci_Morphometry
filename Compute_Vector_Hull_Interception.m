function varargout = Compute_Vector_Hull_Interception(varargin)
%
% Syntax :
%       [intercPoints, order] = Compute_Vector_Hull_Interception(HullSurfMat, P1, P2);
%
% This function computes interception between the vector (P1,P2) and the
% Hull surface.
%
%
% Input Parameters:
%       HullSurfMat             : Hull Surface (Matlab Format)
%       P1                      : First Point of the line segment.
%       P2                      : Second Point of the line segment.
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
if nargin<3 % the indispensable input arguments are not provided
    error('Three inputs are mandatory');
    return
else
    HullSurfMat  = varargin{1};
    P1 = varargin{2};
    P2 = varargin{3};
    HullSurfMat = Surface_Checking(HullSurfMat);
end
%% =========================== Input parameters  =========================%
if nargin > 3
    error('To many inputs');
    return;
end

%% ==================== End of Input parameters  =========================%

%% ================= Computing Interception ============================= %

cont = 0;

% Pial Surface planes
HullSurfMat = Compute_Surface_Normals(HullSurfMat);
VertH = HullSurfMat.SurfData.vertices;FacesH = HullSurfMat.SurfData.faces;NormalsH = HullSurfMat.SurfData.VertexNormals;
np = cross(VertH(FacesH(:,1),:)-VertH(FacesH(:,2),:),VertH(FacesH(:,3),:)-VertH(FacesH(:,2),:),2); Dp = -1*dot(np,VertH(FacesH(:,1),:),2);

Npoints = length(P1);
cont = 0;

tempNorm = normm(P1 - P2);
vectOrient = [P1 - P2]./tempNorm;
vectOrientIndex = [1 2; 2 1];

for po =1:Npoints
    Num = dot(np,repmat(P1(po,:),[size(np,1) 1]),2)+Dp; Den = dot(np,repmat(P2(po,:)-P1(po,:),[size(np,1) 1]),2); % Creating Lines
    
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
    [inte,indord] = unique([xint(ind) yint(ind) zint(ind)],'rows');ind = ind(indord); % Line interception with pial surface
    
    if ~isempty(inte)
        t = NormalsH(FacesH(ind,:),:);
        intNormals = reshape(mean(reshape(reshape(t,[length(ind) prod(size(t))/length(ind)])',[3 length(ind)*3]))',[3 length(ind)])';
        norma = normm(intNormals);
        intNormals = intNormals./[norma norma norma];
        %     orientsign = sign((dot(inte -repmat(Surf2Process.SurfData.vertices(po,:),[size(inte,1) 1]),repmat(Surf2Process.SurfData.VertexNormals(po,:),[size(inte,1) 1]),2))); % Detecting Orientation
        distances = sqrt(sum((inte - repmat(P1(po,:),[size(inte,1) 1])).^2,2)); % Computing distance to the pial surface
        
        [valMin,posit] = min(distances);
        MinPositCoord = inte(posit,:);
        MinPositNorma = intNormals(posit,:);
        
        if sum(MinPositNorma.*vectOrient,2) >0
            [1 2]
        else
            [2 1]
        end
        
    end % End of Looking for non empty values on both sides of the sulci
end
%% ==================== End of Computing Interception =================== %
% Outputs
varargout{1} = Surfwidth;
varargout{2} = sulcalWidth;
if opts.verbose
    toc;
end
return;