function [m1, m2, idx1Out, idx2Out] = match_points_nn(xy1, xy2, maxDist)
%MATCH_POINTS_NN One-way nearest neighbours with unique target assignments.
m1 = zeros(0,2);
m2 = zeros(0,2);
idx1Out = zeros(0,1);
idx2Out = zeros(0,1);
if isempty(xy1) || isempty(xy2)
    return;
end

[idx2Nearest, dist12] = knnsearch(xy2, xy1);
keep = dist12 <= maxDist;
if ~any(keep)
    return;
end

idx1Cand = find(keep);
idx2Cand = idx2Nearest(keep);
distCand = dist12(keep);

fillVal = cast(inf, 'like', distCand);
bestDistByTarget = accumarray(idx2Cand, distCand, [size(xy2,1), 1], @min, fillVal);
isBestForTarget = distCand == bestDistByTarget(idx2Cand);
idx1BestCand = idx1Cand(isBestForTarget);
idx2BestCand = idx2Cand(isBestForTarget);

% Equal-distance ties retain the first source point.
[~, firstPos] = unique(idx2BestCand, 'stable');
idx1Out = idx1BestCand(firstPos);
idx2Out = idx2BestCand(firstPos);

m1 = xy1(idx1Out, :);
m2 = xy2(idx2Out, :);
end
