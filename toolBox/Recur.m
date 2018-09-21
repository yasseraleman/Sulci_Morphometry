function Is = Recur(Is, Surf,level);
%
% Syntax :
% Is = Recur(Is, Surf,level);
%
% This function refills labeled surfaces that contains non-labeled points. It uses
% a recursive process to use neightboor labeled ponts information for non-labeled ones.
%
% Input Parameters:
%   Is       : Surfaces labels.
%   Surf     : Struct variable that contains surface information.
%   level    : Neightborhood level(How many points do we have to take into 
%              account for labelling correction).
%
% Output Parameters:
%  Is        : Corrected Surface labels
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Atlas_Surf
% Plot_oversurf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 30th 2006
% Version $1.0

%=========================Main program====================================%
ind = find(Is == 0);
if ind~=0
    for j = 1: length(ind)
        ind1 =  Surf.Tri(ind(j),3:Surf.Tri(ind(j),2)+2);
        Vert = Surf.SurfData.faces(ind1,:);Vert = unique(Vert(:));ind1 = find(Vert ==ind(j));
        Vert(ind1) = [];
        ord(j) = length(find(Is(Vert) ~=0 ));
    end
    [x,y] = sort(ord,'descend');indtemp = find(x);
    indtemp1 = ind(y(indtemp));clear ind; ind = indtemp1;
    for j = 1:length(ind)
        ind1 =  Surf.Tri(ind(j),3:Surf.Tri(ind(j),2)+2);
        Vert = Surf.SurfData.faces(ind1,:);Vert = unique(Vert(:));ind1 = find(Vert == ind(j));
        Vert(ind1) = [];It = Is(Vert); indt = find(It~=0);
        Vert = [ind(j);Vert(indt)];
        D = dist(Surf.SurfData.vertices(Vert,:)'); D = D(1,2:end);
        if size(D,2)>=level
            indt = find(D == min(D));
            Is(ind(j)) = Is(Vert(indt(1)+1));
        else
            Is(ind(j)) = 0;
        end
    end
    indr = find(Is == 0);
    l = ismember(ind,indr);
    if sum(l(:)) ==size(ind,1)
        return;
    end
    if ~isempty(indr)
        Is = Recur(Is, Surf,level);
    end
end
%========================End of main program==============================%
return;