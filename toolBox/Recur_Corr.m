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
if ~isempty(ind)
    A = unique(Surf.SurfData.faces(Surf.Tri(ind(1),3:2+Surf.Tri(ind(1),2)),:));A = A(:);
    indt = find(Surf.Is(A)~= lab);
    A(indt) = [];A(A == ind(1))= [];
    T = unique([ind(1); A(:)]);
    labid(T) = cont;
    An = rand(size(A));
    while sum(A)~=sum(An)
        An = A;
        Neib = Surf.Tri(A,3:end); Neib = unique(nonzeros(Neib(:)));
        A = unique(Surf.SurfData.faces(Neib,:));A = A(:);
        indt = find(Surf.Is(A)~= lab);
        A(indt) = [];
        labid(A) = cont;
        T =unique([T;A(:)]);
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
            T = [T; A(:)];
        end
    else
        return;
    end
end