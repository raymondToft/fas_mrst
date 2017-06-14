function part = partitionAgglomerationPairwise(A, coupleTol, checkDD, varargin)
    assert((coupleTol <= 1 && coupleTol >= 0) || isinf(coupleTol))
    % Transform to positive, remove diagonal
    D = diag(A);
    
    A_offdiag = abs(A - diag(D));
    connStr = full(max(A_offdiag, [], 2));
    
    % Undecided variables
    U = true(size(A,1),1);
    if checkDD
        % check for diagonal dominance
        U(5*sum(A_offdiag, 2) < D) = false;
    end
    
    [rowi, coli, conn] = find(A_offdiag);
    
    if numel(varargin)
        % Strict tolerance
        isStrong = conn >= varargin{1};
    else
        isStrong = conn >= coupleTol*connStr(rowi);
    end
    
    measure =  accumarray(rowi, isStrong);
    n_c = 1;
    C = nan(size(A, 1), 2);
    
    while any(U)
        [v, nextCoarse] = min(measure); %#ok
        targetConn = rowi == nextCoarse;
        candidates = coli(targetConn);
        lconn      = conn(targetConn); 
        [v, pair_loc] = max(lconn); %#ok
        pair = candidates(pair_loc);
        
        c = nextCoarse;
        if isStrong(targetConn & coli == pair)
            c = [c, pair]; %#ok not really the case.
        end
        
        locstrong = [];
        for i = 1:numel(c)
            C(n_c, i) = c(i);
            locstrong = [locstrong; getStrong(c(i), coli, rowi, isStrong)]; %#ok
            conn(coli == c(i)) = 0;
            isStrong(coli == c(i)) = 0;
        end
        measure(locstrong) = max(measure(locstrong) - 1, 0);
        measure(c) = nan;
        measure(measure == 0) = nan;
        U(c) = false;
        
        n_c = n_c + 1;
        
        if all(isnan(measure))
            % We might wind up with the odd man out, assign explicitly
            if any(U)
                last = find(U);
                assert(numel(last) <= 2);
                C(n_c, 1:numel(last)) = last;
            end
            break
        end
        if mod(n_c, 100) == 0
            dispif(mrstVerbose, '%d of maximally %d\n', n_c, size(A, 1));
        end
    end
    
    part = zeros(size(A, 1), 1);
    for i = 1:2
        v = ~isnan(C(:, i));
        part(C(v, i)) = find(v);
    end
end

function s = getStrong(i, coli, rowi, isStrong)
    s = coli(isStrong(rowi == i));
end