function varargout = Sort_3D_pointCloud(varargin);
%
% Syntax :
%       [maxdistPair] = Sort_3D_pointCloud(spacCoord);
%
% This function computes point cloud extremes according to a distance
% function.
%
%
% Input Parameters:
%       spacCoord               : 3D coordinates
%
% Output Parameters:
%       maxdistPair             : Pair of points with high distance between
%                                 them
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% February 13th 2015
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin < 1
    error('Please enter points coordinates');
    return
end
if nargin == 1
    spacCoord = varargin{1};
    method = 'mixed';
end
if nargin == 2
    spacCoord = varargin{1};
    method = varargin{2};
    switch   lower(deblank(method))
        case 'euclidean'
        case 'geodesic'
        case 'mixed'
        otherwise
            error('Only euclidean or geodesic can be selected as methods');
            return;
    end
end
%% ====================== End of Input parameters  =======================%




%% =========================== Main Program ============================= %
if size(spacCoord,2)>3
    clustId = unique(nonzeros(spacCoord(:,4)));
else
    clustId = 1;
    spacCoord(:,4) = ones(size(spacCoord,1),1);
end

Nc = length(unique(clustId));
maxdistPair = [0 0 0];
for i = 1:Nc
    ind = find(spacCoord(:,4) == clustId(i));
    spacCoordCl = spacCoord(ind,:);
    % Creating Auxiliar order Matrix
    Npoints = size(spacCoordCl,1);
    pointMat = repmat([1:Npoints],[Npoints 1])';
    pointMat(logical(eye(size(pointMat)))) = [];
    
    ordMat = reshape(pointMat,[Npoints-1 Npoints])'; % Order Matrix for processing
    
    % Points cloud coordinates
    Xinterp = spacCoordCl(:,1);
    Yinterp = spacCoordCl(:,2);
    Zinterp = spacCoordCl(:,3);
    
    % Differences Matrix
    Xdiff = Xinterp(ordMat) - repmat(Xinterp([1:Npoints]'),[1 Npoints-1]);
    Ydiff = Yinterp(ordMat) - repmat(Yinterp([1:Npoints]'),[1 Npoints-1]);
    Zdiff = Zinterp(ordMat) - repmat(Zinterp([1:Npoints]'),[1 Npoints-1]);
    
    % Methods to find point cloud extrems
    switch method
        case 'euclidean'
            [~,loc] = max(sqrt(Xdiff(:).^2 + Ydiff(:).^2 + Zdiff(:).^2));
            [xi,yi] = ind2sub(size(ordMat),loc);
            maxdistPair= [maxdistPair;ind(xi),ind(ordMat(xi,yi)) clustId(i)];
        case 'geodesic'
            eucDist = dist([Xinterp Yinterp Zinterp]');
            geoDist= sparse(length(Xinterp),length(Xinterp));
            index = sub2ind(size(geoDist) ,[1:length(Xinterp)]' ,[2:length(Xinterp) 1]');
            geoDist(index) = eucDist(index);
            geoDist = geoDist + geoDist';
            geoDist = graphallshortestpaths(geoDist);
            [~,loc] = max(geoDist(:));
            [xi,yi] = ind2sub(size(geoDist),loc);
            maxdistPair= [maxdistPair;ind(xi),ind(yi) clustId(i)];
        case 'mixed'
            eucDist = dist([Xinterp Yinterp Zinterp]');
            geoDist= sparse(length(Xinterp),length(Xinterp));
            index = sub2ind(size(geoDist) ,[1:length(Xinterp)]' ,[2:length(Xinterp) 1]');
            geoDist(index) = eucDist(index);
            geoDist = geoDist + geoDist';
            geoDist = graphallshortestpaths(geoDist);
            [~,loc] = max(geoDist(:)+eucDist(:));
            [xi,yi] = ind2sub(size(geoDist),loc);
            maxdistPair= [maxdistPair;ind(xi),ind(yi) clustId(i)];
    end
%     [xi,yi] = ind2sub(size(ordMat),loc);
%     maxdistPair= [maxdistPair;ind(xi),ind(ordMat(xi,yi)) clustId(i)];
end
maxdistPair(1,:) = [];

%% ==================== End of Main Program ============================= %
varargout{1} = maxdistPair;

return;