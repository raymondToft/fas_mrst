function [restrict, prolong] = tensorLinearInterpolator(A, q, p, isFW, varargin)
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

    is3D = G.griddim == 3;
    
    if ~is3D
        [ii, jj] = gridLogicalIndices(G);
        kk = ones(G.cells.num, 1);
    else
        [ii, jj, kk] = gridLogicalIndices(G);
    end
    % Renum so ijk is in column first indexing scheme
    ijk = {jj, ii, kk};
    
    tmp = cell(CG.faces.num, 1);
    
    
    % Begin by interpolating over all coarse faces, giving the linear
    % functions from each "center" to next center
    faceInterp = struct('direction', tmp, 'values', tmp, 'cells', tmp);
    for i = 1:CG.faces.num
        N = CG.faces.neighbors(i, :);
        
        dirmap = [2 1 3];
        if any(N == 0)
            faceInterp(i).direction = -1;
            [a, b] = deal(p.centers(N(N~=0)));
            
            if N(1)
                x1 = CG.cells.centroids(N(1), :);
                x2 = CG.faces.centroids(i, :);
            else
                x1 = CG.faces.centroids(i, :);
                x2 = CG.cells.centroids(N(2), :);
            end
            tmp = x2 - x1;
            [~, iii] = max(abs(tmp));
            direction = dirmap(iii);
            if sign(CG.faces.normals(i, iii)) < 0
                % Switch points depending on normal component of the face
                tmp = x1;
                x1 = x2;
                x2 = tmp;
            end
        else
            centers = p.centers(N);
            a = centers(1); 
            b = centers(2);
            
            x1 = CG.cells.centroids(N(1), :);
            x2 = CG.cells.centroids(N(2), :);
            if ii(a) ~= ii(b)
                direction = 2;
            elseif jj(a) ~= jj(b)
                direction = 1;
            else
                direction = 3;
            end
        end
        
        cells = true(G.cells.num, 1);
        for j = 1:G.griddim
            if j == direction
                % Include centers in subset
                d = dirmap(direction);
                c = G.cells.centroids(:, d);
                cells = cells & c >= x1(d) & c <= x2(d);
            else
                cells = cells & (ijk{j} == ijk{j}(a) | ijk{j} == ijk{j}(b));
            end
        end

        
        cells = find(cells);
        
        x = G.cells.centroids(cells, :);

        v = interpolate(x1, x2, x);

        s = struct('direction', direction, ...
                   'values',    v,         ...
                   'cells',     cells);
       
        faceInterp(i) = s;

    end
    
    [v_i, v_j, v_c] = deal([]);
    
    for i = 1:CG.cells.num
        faces = gridCellFaces(CG, i);
        
        facei = faceInterp(faces);
        faceDir = arrayfun(@(x) x.direction, facei);
        
        nDir = max(faceDir);
        
        values = cell(nDir, 1);
        indexSpace = cell(nDir, 1);
        for j = 1:nDir
            isjdir = find(faceDir == j);
            if numel(isjdir) > 2
                normals = CG.faces.normals(faces(isjdir), :);
                [~, subs] = unique(normals, 'rows');
                isjdir = isjdir(subs);
            end
            fi = facei(isjdir);
            fa = faces(isjdir);
            [v, cells] = getCellInterp(G, CG, i, fi, fa, j);
            if numel(v) == 1 && v == 0
                % Avoid collapse in 1D...
                v = 1;
            end
            values{j} = v;

            indexSpace{j} = ijk{j}(cells);
        end
        % Use tensor product to define the values outside the linear
        % segments, then do some reshape magic to get index space blocks
        v = tensorProduct(values);
        n = numel(G.cartDims);
        for k = 1:n
            tmp = size(v);
            if numel(tmp) == 2;
                tmp = [tmp, 1]; %#ok 
            end
            tmp(k) = 1;

            dims = {1, 1, 1};
            dims{k} = [];
            if n < 3
                d = reshape(indexSpace{k}, dims{1}, dims{2});
                indexSpace{k} = repmat(d, tmp(1), tmp(2));
            else
                d = reshape(indexSpace{k}, dims{1}, dims{2}, dims{3});
                indexSpace{k} = repmat(d, tmp(1), tmp(2), tmp(3));
            end
            
        end
        subs = sub2ind(G.cartDims, indexSpace{dirmap(1:numel(G.cartDims))});
        
        % A cell is inside if present in the indexmap list, else virtual
        inside = ismember(subs, G.cells.indexMap);
        
        % Flatten
        subs = subs(inside(:));
        v = v(inside(:));
        
        % Remove any items represented more than once
        [subs, si] = unique(subs);
        v = v(si);
        
        % Enforce v in (0, 1)
        v = max(v, 0);
        v = min(v, 1);
        
        cells = (find(ismember(G.cells.indexMap, subs)));
        
        v_i = [v_i; cells]; %#ok
        v_j = [v_j; repmat(i, numel(cells), 1)]; %#ok
        v_c = [v_c; v(:)]; %#ok
    end

    n = numel(partition);
    m = max(partition);
    
    prolong = sparse(v_i, v_j, v_c, n, m);
    if isFW
        % Ensure that interpolation is ok even for e.g. 1 cell
        sump = sum(prolong, 2);
        bad = sump == 0;
        sump(bad) = .1;
        
        prolong = bsxfun(@rdivide, prolong, sump);
        prolong = prolong + sparse(find(bad), CG.partition(bad),...
                        ones(sum(bad), 1), G.cells.num, CG.cells.num);
    end
    restrict = prolong .';
end

function v = interpolate(p1, p2, pv)
    v = sqrt(sum( bsxfun(@minus, p2, pv).^2, 2))./norm(p2-p1, 2);
end

function [v, cells] = getCellInterp(G, CG, cellIndex, fi, faces, dir)
    tmp = nan(G.cells.num, 1);
    for i = 1:numel(fi)
        first = CG.faces.neighbors(faces(i), 1);
        if (first == cellIndex) || (first ~= cellIndex && dir == 3)
            tmp(fi(i).cells) = fi(i).values;
        else
            tmp(fi(i).cells) = 1 - fi(i).values;
        end
    end
    ok = ~isnan(tmp);
    
    v = tmp(ok);
    cells = find(ok);
%     cells = G.cells.indexMap(find(ok));
end

function v = tensorProduct(values)
    nv = numel(values);
    
    v = values{1};
    if nv > 1
        v = kron(v, values{2}');
    end

    if nv > 2
        % Fallback because not all MATLAB versions support n dimensional
        % kronecker product
        N = numel(values{3});

        tmp = zeros(size(v, 1), size(v, 2), N);
        for i = 1:size(v, 1)
            for j = 1:size(v, 2)
                for k = 1:N
                    tmp(i, j, k) = v(i, j)*values{3}(k);
                end
            end
        end
        v = tmp;
    end
end
