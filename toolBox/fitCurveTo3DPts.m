function [cPts] = fitCurveTo3DPts(XYZ, stPt, endPt, Npoints, iter);
MAX_ITERATIONS = 5;
completedIterations = 0;

% Initialise control pts linearly between start/end anchors
cPts = interp1(0:1, [stPt; endPt], linspace(0,1,Npoints));

while completedIterations < MAX_ITERATIONS
    % Get the nearest-neighbour cntrl pt for each of the sample points
    sqDists = cellfun(@(x)sum(bsxfun(@minus, XYZ, x).^2,2), num2cell(cPts,2),'UniformOutput', false);
    [~, nnIdxs] = min(cat(2,sqDists{:}),[],2);
    % Keep the anchors, update inner cPts to the mean of their linked input pts
    for i = 2:Npoints-1
        cPts(i,:) = mean(XYZ(nnIdxs==i,:),1);
    end
    % Handle any cPts that didn't have linked pts so their mean became NaN
    goodIdxs = find(~isnan(cPts(:,1)));
    badIdxs = find(isnan(cPts(:,1)));
    cPts(badIdxs,:) = interp1(goodIdxs, cPts(goodIdxs,:), badIdxs);
    % Re-spread the control points out linearly
    cPtCumSumDists = cumsum([0; sqrt(sum(diff(cPts,1).^2,2))]);
    cPts = interp1(cPtCumSumDists, cPts, linspace(0, cPtCumSumDists(end), Npoints));
    completedIterations = completedIterations + 1;
end


% Smoothing
P1 = [cPts(1:size(cPts,1)-2,:)];
P2 = [cPts(2:size(cPts,1)-1,:)];
P3 = [cPts(3:size(cPts,1),:)];
n = cross(P1-P2,P3-P2,2); D = -1*dot(n,P1,2);
t = P2-P1; u = P3-P1; v = P3-P2;
t2 = sum((t.^2)')'; u2 = sum((u.^2)')'; w2 = sum((n.^2)')';
c = P1+(repmat(t2.*dot(u,v,2),[1 size(u,2)]).*u-repmat(u2.*dot(t,v,2),[1 size(t,2)]).*t)./(2*repmat(w2,[1 size(P1,2)])+eps);
r = 1/2*sqrt(t2.*u2.*dot(v,v,2)./w2+eps);
curv = ones(length(r),1)./(r + eps);
ind = find(curv > 4);
ind = ind +1;
endpoint = cPts(end,:);
initpoint = cPts(1,:);
cPts(ind,:) = [];
%iter = 0;
% if ~isempty(ind)
%     if (iter > 5)
%         return;
%     end
    while (iter <= 10)&~isempty(ind)
        iter  = iter + 1;
        cPts = fitCurveTo3DPts(cPts, initpoint, endpoint, Npoints, iter);
        P1 = [cPts(1:size(cPts,1)-2,:)];
        P2 = [cPts(2:size(cPts,1)-1,:)];
        P3 = [cPts(3:size(cPts,1),:)];
        n = cross(P1-P2,P3-P2,2); D = -1*dot(n,P1,2);
        t = P2-P1; u = P3-P1; v = P3-P2;
        t2 = sum((t.^2)')'; u2 = sum((u.^2)')'; w2 = sum((n.^2)')';
        c = P1+(repmat(t2.*dot(u,v,2),[1 size(u,2)]).*u-repmat(u2.*dot(t,v,2),[1 size(t,2)]).*t)./(2*repmat(w2,[1 size(P1,2)])+eps);
        r = 1/2*sqrt(t2.*u2.*dot(v,v,2)./w2+eps);
        curv = ones(length(r),1)./(r + eps);
        ind = find(curv > 4);
        ind = ind +1;
        cPts(ind,:) = [];
    end
% end

% end of smoothing
return;