function [restrict, prolong] = linearInterpolator(A, q, p, varargin)
    cache = true;
    if numel(varargin)
        cache = varargin{1};
    end
    
    partition = p.partition;
    CG = p.grid;
    G = CG.parent;
    
    if isfield(p, 'R') && isfield(p, 'I') && cache
        %caching for speedup
        restrict = p.R;
        prolong = p.I;
        return
    end
    

    
    N = CG.faces.neighbors;
    
    coarseCellNum = accumarray(partition, ones(G.cells.num, 1));
    coarseVals = arrayfun(@(x) nan(x, 1), coarseCellNum, 'UniformOutput', false);
    [v_i, v_j, v_c] = deal([]);
    
    ijk = gridLogicalIndices(CG);
    
    for i = 1:CG.cells.num
        fa = gridCellFaces(CG, i);

        if isfield(G, 'cartDims') && 0
            
%             isLocalCoarse = false(CG.cells.num, 1);
            isLocalCoarse = false(CG.cells.num, 1);
            for ii = -1:1
                for jj = -1:1
                    if G.griddim == 3
                        for kk = -1:1
                            isLocalCoarse = isLocalCoarse | ...
                                            (ijk{1} == ijk{1}(i) + ii & ...
                                             ijk{2} == ijk{2}(i) + jj & ...
                                             ijk{3} == ijk{3}(i) + kk);
                        end
                    else
                        isLocalCoarse = isLocalCoarse | ...
                                        (ijk{1} == ijk{1}(i) + ii & ...
                                         ijk{2} == ijk{2}(i) + jj);
                    end
                end
            end
            
            N = find(isLocalCoarse);
            N = N(N~=i);
        else
            N = CG.faces.neighbors(fa, :);
            isInner = all(N~=0, 2);
            N = N(isInner, :);        
            N = unique(N(:));
            N = N(N~=i);
        end
%         m = commonNeighbors(CG, N);
%         N = [N; m];
%         N = unique(N);

        % Put current node first
        N = [i; N];
        
        
        cells = find(ismember(partition, unique(N(:))));
        cent = G.cells.centroids(cells, :);
        values = zeros(size(cells));


        
        dpt = zeros(size(N));
        dpt(1) = 1;
        
        c = G.cells.centroids(p.centers(N), :);
        
        if numel(dpt) == 1
            values = ones(size(cells));
        else
            F = scatteredInterpolant(c, dpt);
            values = F(cent);
            
            values(values < 0) = 0;
            values(values > 1) = 1;
        end
%         if numel(N) > 0
%             figure(1); clf;
%             plotCellData(G, values, cells);
%             outlineCoarseGrid(G, CG.partition)
%             plotGrid(G, cells, 'facec', 'none', 'edgec', 'w', 'edgea', .2);
%             colorbar
%             plotGrid(G, p.centers(N), 'facec', 'r', 'facea', .5)
%             pause;
%         end
        v_i = [v_i; cells]; %#ok
        v_j = [v_j; repmat(i, numel(cells), 1)]; %#ok
        v_c = [v_c; values]; %#ok
    end

    n = numel(partition);
    m = max(partition);
    
    prolong = sparse(v_i, v_j, v_c, n, m);
%     prolong = bsxfun(@rdivide, prolong, sum(prolong, 2));
    restrict = prolong .';
end

function m = commonNeighbors(G, n)
    N = G.faces.neighbors;
    m = G.faces.neighbors(ismember(N(:, 1), n) & ismember(N(:, 2), n), :);
    m = unique(m(:));
end