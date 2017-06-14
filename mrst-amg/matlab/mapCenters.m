function centers = mapCenters(CG, G)
    centers = zeros(CG.cells.num, 1);
    for i = 1:CG.cells.num
        local = find(CG.partition == i);
        
        dist = bsxfun(@(x, y) x-y, G.cells.centroids(local,:), ...
                                   CG.cells.centroids(i,:));
        [v, ind] = min(sqrt(sum(dist.^2, 2))); %#ok
        centers(i) = local(ind);
    end
end