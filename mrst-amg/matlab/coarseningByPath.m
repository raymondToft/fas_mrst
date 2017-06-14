function partition = coarseningByPath(A, G, partition, L, P)
    Nc = sum(partition.coarse);

    p = sparse(find(partition.coarse), 1:Nc, ones(Nc, 1), size(A, 1), Nc);

    N = getNeighbourship(G);
    A_n = getConnectivityMatrix(N);
    % Use matrix multiplication to efficiently obtain distance between
    % nodes
    for i = 1:L
        p = p + A_n * p;
    end
    
    % Define connections between coarse blocks via paths, and select the
    % strong ones via the P threshold.
    p = p(partition.coarse, :);
    p = p .';
    coarse_partition = standardCoarsening(p, inf, P);
    
    % Remap to original ordering
    oldcoarse = find(partition.coarse);
    newcoarse = false(G.cells.num, 1);
    newcoarse(oldcoarse(coarse_partition.coarse)) = true;
    newfine = ~newcoarse;
    
    % Store and return
    partition.coarse = newcoarse;
    partition.fine = newfine;
    
    
    nc = numel(coarse_partition.strongCoarseConnected);
    newconn = coarse_partition.strongCoarseConnected;
    for i = 1:nc
        self = coarse_partition.strongCoarseConnected{i};
        old = partition.strongCoarseConnected(self);
        newconn{i} = vertcat(old{:});
    end
    partition.strongCoarseConnected = newconn;
    
end