function partitions = makePartitions(G, A, q, getOperators, coarsening, addTracer, opt)
    % Set up initial dimensions to be equivialent to the fine scale
    dims = G.cartDims;
    
    if numel(coarsening) == 1
        coarsening = repmat(coarsening, 1, numel(dims));
    end
    
    partitions = {};
    converged = false;
    iter = 1;
    while ~converged && iter <= opt.maxCoarseLevels
        % Ensure that dimensions are always more than zero and that the
        % grid is not coarsened to a finer degree than the grid dimensions
        dims = min(max(ceil(dims./coarsening), 1), G.cartDims);
        
        switch lower(opt.coarsening)
            case 'uniform'
                p = partitionUI(G, dims);
                converged = ~any(dims > 2);
            case 'metis'
                if prod(dims) == 1
                    p = ones(G.cells.num, 1);
                else
                    p = callMetisMatrix(A, prod(dims), ...
                        '-minconn -contig -ufactor=500  -no2hop -objtype=cut');
                end
                % Maybe ok??
                converged = ~any(dims > 2);
            case 'pairwise'
                p = partitionAgglomerationLevels(A, opt.cfstrength, iter == 1, opt.levels);
                converged = max(p) <= max(opt.maxCoarseSize, 2) || 1.25*max(p) > size(A, 1);
            otherwise
                error('Unknown coarsening strategy');
        end
        p = processPartition(G, p);
        p = compressPartition(p);
        G_old = G;
        
        G = generateCoarseGrid(G, p);

        % Add fields to the coarse grid to ensure that it passes as a
        % regular grid for our purposes.
        G.cartDims = dims;
        
        
        G = coarsenGeometry(G);
        centers = mapCenters(G, G.parent);
        
%          oldmap = G_old.cells.indexMap(centers);
%          ijk = gridLogicalIndices(G_old, oldmap);
         
         G.cells.indexMap = remapLogicalIndices(G_old.cartDims, G.cartDims,...
                                                G_old.cells.indexMap(centers), coarsening);
        clear G_old
        
%         G.cells.indexMap = (1:G.cells.num) .';
        if addTracer
            G = addTracerData(G);
        end
        
        if isfield(G.parent, 'parent')
            G.parent = rmfield(G.parent, 'parent');
        end
        
        partitions = [partitions; ...
            struct('partition', p, 'grid', G, 'centers', centers)]; %#ok
        
        % Get partitions and store them in the array
        [R, I] = getOperators(A, q, partitions{end});
        if opt.jacobiBasis
            I = jacobiUpdateBasis(A, I, opt.jacobiBasis, 2/3, opt.cfstrength);
            I = bsxfun(@rdivide, I, sum(I, 2));
            R = I';
        end
        partitions{end}.R = R;
        partitions{end}.I = I;
        
        % Store the reduced systems
        A = R*A*I; q = R*q;
        
        partitions{end}.A = A;
        partitions{end}.q = q;
        
        dispif(opt.verbose, 'Partitioned level %d of maximum %d\n', iter, opt.maxCoarseLevels);

        iter = iter + 1;
    end
end

function G = addTracerData(G)
    % Add dummy fields to allow for mstpfa-like multigrid method
    G.faces.nodePos = G.parent.faces.nodePos;
    G.faces.nodes = G.parent.faces.nodes;
    G.nodes.coords = G.parent.nodes.coords;
    G.nodes.num = G.parent.nodes.num;
    G = partitionCoarseFaceRegions(G);
    G = partitionMultiscaleTracer(G, [], 'UniformRock', true);
end

function newindex = remapLogicalIndices(olddims, newdims, oldindex, coarsening)
    % Map the indices to the (i,j,k) space, accounting for variable
    % dimensionality.
    [oldSubs{1:numel(olddims)}] = ind2sub(olddims, double(oldindex));
    newSubs = oldSubs;
    % For each ijk direction, coarsen it based on the coarsening factors
    for i = 1:numel(oldSubs)
        newSubs{i} = min(max(ceil(oldSubs{i}./coarsening(i)), 1), olddims(i));
    end
    % Map back to linear indexes of some Cartesian block
    newindex = sub2ind(newdims, newSubs{:});
end
