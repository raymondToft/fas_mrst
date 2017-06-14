function [restrict, prolong] = simpleLinearInterpolator(A, q, p, isFW, varargin)
    cache = true;
    if numel(varargin)
        cache = varargin{1};
    end
    
    partition = p.partition;
    CG = p.grid;
    G = CG.parent;
    
    assert(all(mod(G.cartDims, 2)) == 0)
    assert(all(mod(G.cartDims./CG.cartDims, 2)) == 0)
    
    
    if isfield(p, 'R') && isfield(p, 'I') && cache
        %caching for speedup
        restrict = p.R;
        prolong = p.I;
        return
    end
    
    [v_i, v_j, v_c] = deal([]);
    
    is3D = G.griddim == 3;
    
    if ~is3D
        [ii, jj] = gridLogicalIndices(G);
        kk = ones(G.cells.num, 1);
    else
        [ii, jj, kk] = gridLogicalIndices(G);
    end
    
    for i = 1:CG.cells.num
        values = getValues(is3D, isFW);
        tmp = repmat([-1 0 1], 3, 1);
        
        i_ind = zeros(size(values));
        for j = 1:size(values, 3)
            i_ind(:, :, j) = tmp;
        end
%         i_ind = repmat(tmp, 1, 1, 1 + 2*(G.griddim == 3));
        
        if ~is3D            
            j_ind = i_ind .';
            k_ind = zeros(3, 3);
        else            
            j_ind = permute(i_ind, [1 3 2]);
            k_ind = permute(i_ind, [2 3 1]);
        end
        ipos = ii(p.centers(i)) - i_ind;
        jpos = jj(p.centers(i)) - j_ind;
        kpos = kk(p.centers(i)) - k_ind;

        ipos = ipos(:);
        jpos = jpos(:);
        kpos = kpos(:);

        values = values(:);

        nc = numel(ipos);

        cells = nan(nc, 1);

        for j = 1:nc
            if ipos(j) && jpos(j) && kpos(j);
                cells(j) = find(ii == ipos(j) & ...
                                jj == jpos(j) & ...
                                kk == kpos(j));
            end
        end
        inside = ~isnan(cells);
        cells = cells(inside);
        values = values(inside);


        assert(numel(cells) == numel(unique(cells)))

        v_i = [v_i; cells]; %#ok
        v_j = [v_j; repmat(i, numel(cells), 1)]; %#ok
        v_c = [v_c; values]; %#ok
    end

    n = numel(partition);
    m = max(partition);
    
    prolong = sparse(v_i, v_j, v_c, n, m);
    if isFW
        % Handle boundary brute force
        prolong = bsxfun(@rdivide, prolong, sum(prolong, 2));
    end
    restrict = prolong .';
end

function values = getValues(is3D, isFW)
    if is3D
        values = zeros(3, 3, 3);
        if isFW
            edge = [1, 2, 1; ...
                    2, 4, 2; ...
                    1, 2, 1];

            middle = [2, 4, 2; ...
                      4, 8, 4; ...
                      2, 4, 2];
        else
            edge = [0, 0, 0; ...
                    0, 1, 0; ...
                    0, 0, 0];
            
            middle = [0, 1, 0; ...
                      1, 6, 1; ...
                      0, 1, 0];
        end
        values(:, :, 1) = edge;
        values(:, :, 2) = middle;
        values(:, :, 3) = edge;
    else
        if isFW
            values = [1, 2, 1; ...
                      2, 4, 2;...
                      1, 2, 1];
        else
            values = [0, 1, 0; ...
                      1, 4, 1;...
                      0, 1, 0];
        end
    end
    values = values./sum(values(:));
end