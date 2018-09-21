function [Surf] = Surf_Corr(Surf);
%
% Syntax :
% [Surf] = Surf_Corr(Surf);
%
% This function corrects labeled surfaces. It removes 
% small isolated points( or group of points) that are 
% surrounded by a different label.     
%
% Input Parameters:
%   Surf      : Individual Surfaces.
%
% Output Parameters:
%   Surf  : Corrected Surfaces.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Plot_oversurf 
% Exp_Surf
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1st 2007
% Version $1.0

% [fa,intt] =sort(double(Surf.SurfData.faces)');
% [ford,int2] =sortrows(fa');
% temp = ford(1:end-1,:) - ford(2:end,:);
% fnd = find((temp(:,1)== 0)&(temp(:,2)== 0)&(temp(:,3)== 0));
% Surf.SurfData.faces(int2(fnd),:)=[];
% Surf.SurfData.faces =uint32(Surf.SurfData.faces);Nfaces =size(Surf.SurfData.faces,1);
% [Tri] = Vert_Neib(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),Nfaces);Surf.Tri=Tri;
% Temp = sum(Tri);
% Tri(:,Temp==0) = [];
strl = unique(Surf.Is);
if ~isfield(Surf,'Tri')
    Npoints = size(Surf.SurfData.vertices,1);
    Nfaces = size(Surf.SurfData.faces,1);
    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
    Temp = sum(Tri);
    Tri(:,Temp==0) = [];
    Surf.Tri = Tri; clear Tri;
end
lab = Surf.Is;
set(0,'RecursionLimit',1000);
for i = 1:size(strl,1)
    lab = strl(i);
    [labid] = Recur_Corr(Surf,lab,zeros(size(Surf.Is)),1);
    c = accumarray(nonzeros(labid),ones(length(nonzeros(labid)),1));
    indpos = find(c == max(c));
    if length(indpos)==1
        Surf.Is((labid~= indpos(1))&(labid~=0)) = 0;
    else
        Surf.Is(~(ismember(labid,indpos))&(labid~=0)) = 0;
    end
end
Is = Recur(Surf.Is, Surf,1);
Surf.Is = Is;
ind = find(Surf.Is ==0);
if ~isempty(ind)
    dfac = unique(nonzeros(Surf.Tri(ind,3:end)));
    Surf.SurfData.faces(dfac,:) = [];
    Surf.SurfData.vertices(ind,:) = [];
    if isfield(Surf.SurfData,'VertexNormals')
        Surf.SurfData.VertexNormals(ind,:) = [];
    end
    if isfield(Surf,'Is')
        Surf.Is(ind,:) = [];
    end
    if isfield(Surf.SurfData,'FaceVertexCData')
        Surf.SurfData.FaceVertexCData(ind,:) = [];
    end
    cont = 0;
    Mat = Surf.SurfData.faces;
    for i =1:size(ind,1)
        cont= cont+1;
        dvert = find(Surf.SurfData.faces >ind(i));
        Mat(dvert) = Surf.SurfData.faces(dvert) - cont;
    end
    Surf.SurfData.faces = Mat;
    Npoints = size(Surf.SurfData.vertices,1);
    Nfaces = size(Surf.SurfData.faces,1);
    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
    Temp = sum(Tri);
    Tri(:,Temp==0) = [];
    Surf.Tri = Tri;
end



function [labid] = Recur_Corr(Surf,lab,labid,cont);
%
% Syntax :
% [labid] = Recur_Corr(Surf,lab,labid,cont);
%
% Recursive function for labelling correction.    
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Plot_oversurf 
% Exp_Surf
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1st 2007
% Version $1.0

ind = find(Surf.Is ==lab&(labid ==0));
A = unique(Surf.SurfData.faces(Surf.Tri(ind(1),3:2+Surf.Tri(ind(1),2)),:));
indt = find(Surf.Is(A)~= lab);
A(indt) = [];A(A == ind(1))= [];
T = unique([ind(1); A]);
labid(T) = cont;
An = rand(size(A));
while sum(A)~=sum(An)
    An = A;
    Neib = Surf.Tri(A,3:end); Neib = unique(nonzeros(Neib(:)));
    A = unique(Surf.SurfData.faces(Neib,:));
    indt = find(Surf.Is(A)~= lab);
    A(indt) = [];
    labid(A) = cont;
    T =unique([T;A]);
end
indn = find(Surf.Is ==lab&(labid ==0));
if ~isempty(indn)
    cont = cont+1;
    [labid] = Recur_Corr(Surf,lab,labid,cont);
else
    return;
end
if ~isempty(A)
    for i = 1:size(A,2)
        TA = unique(Surf.SurfData.faces(Surf.Tri(A(i),3:2+Surf.Tri(A(i),2)),:));
        indt = find(Surf.Is(TA(i))~= lab);
        TA(indt) = [] ;
        T = [T; A];
    end
else
    return;
end