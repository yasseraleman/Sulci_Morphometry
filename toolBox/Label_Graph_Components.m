function varargout = Label_Graph_Components(varargin);
%
% Syntax :
%     varargout = Label_Graph_Components(varargin);
%
% This script labels the unconnected networks in a connectivity matrix 
% provided they are higher than a certain threshold
%
% Input Parameters:
%     Mat                      : Connectivity Matrix
%     thresh                   : Connectivity Threshold 
%
% Output Parameters:
%   LabNet                     : Labeled connectivity Matrix.
%   modClust                   : Regions that belong to each cluster 
%
%__________________________________________________
% Authors:  Yasser Aleman Gomez 
% LIM
% November 15th 2014
% Version $1.0

%% ===================== Checking Input Parameters ===================== %%
if nargin == 1
    Mat = varargin{1};
end
if nargin == 2
    Mat = varargin{1};
    thresh = varargin{2};
    ind = find(abs(Mat)< thresh);
    Mat(ind) = 0;
end
if nargin > 2
    error('Too Many Inputs');
    return;
end
if nargout > 2
    error('Too Many Outputs');
    return;
end

%% ================== End of Checking Input Parameters ================= %%

% --------------- Example 
% Mat = zeros(20,20);
% Mat(10,1) =1;
% Mat(1,10) =1;
% Mat(10,5) =1;
% Mat(5,10) =1;
% Mat(17,18) =1;
% Mat(18,17) =1;
% Mat(5,4) =1;
% Mat(4,5) =1;

%% ========================== Main Program ============================= %%
% ---------------------- Making a symetric matrix 
Mat(find(speye(size(Mat)))) = 0;
Mat = logical(Mat + Mat');

% -------------- Looking for nonzero connections
ind = find(Mat);
[X,Y] = ind2sub(size(Mat),ind);
LabNet = sparse(size(Mat,1),size(Mat,2));

% ---------- Detecting Graph Components 
[Nclust,allclust] = graphconncomp(sparse(logical(Mat)));
allclust = allclust(:);
temp = accumarray(allclust,allclust*0+1);
clust = find(temp>=2);
Nclust = length(clust);

% ----------- Sorting cluster according to their size 
[tempSort,order] = sort(temp,'descend');
tempSort = tempSort(1:Nclust);
order = order(1:Nclust);

% ----------- Reordering cluster labels
[temp,modClust] = ismember(allclust,order);

% --------- Assigning new cluster Ids
LabNet(ind) = modClust(Y); 

%%  ------------------ Output parameters -------------------------------- %
varargout{1} = LabNet;   % Labelled Connectivity Matrix
varargout{2} = modClust; % Regional Ids

return;
