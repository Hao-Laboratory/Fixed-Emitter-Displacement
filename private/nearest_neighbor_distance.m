function nnDist = nearest_neighbor_distance(xy)
N = size(xy,1);
if N < 2
    nnDist = inf(N,1);
    return;
end

[~, D] = knnsearch(xy, xy, 'K', 2);
nnDist = D(:,2);
end
