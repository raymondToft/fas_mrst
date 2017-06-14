function [restrict, prolong] = multigridInterpolator(A, q, part, useStandard, varargin)
    cache = true;
    if numel(varargin)
        cache = varargin{1};
    end
    if isfield(part, 'R') && isfield(part, 'I') && cache
        %caching for speedup
        restrict = part.R;
        prolong = part.I;
        return
    end
    
    % Get diagonal and connection strengths, ensure all connections to be
    % positive.
    D = diag(A);
    A_offdiag = abs(A - diag(D));
    [rowi, coli, conn] = find(A_offdiag);
    
    % Get the partition indices
    fine   = find(part.fine);
    coarse = find(part.coarse);
    
    % Interpolation of coarse values should always be 1, add to sparse
    % index vectors.
    i_ind = coarse;
    j_ind = (1:numel(coarse)) .';
    v     = ones(numel(coarse), 1);

    % Set up function handle for getting connections, allowing for ignoring
    % of specific connections
    getConn = @(f, ignore) getConnections(f, ignore, rowi, coli, conn, coarse, part);
    for i = 1:numel(fine)
        f = fine(i);
        [N, P, N_conn, P_conn, coarseIndex, local] = getConn(f, nan);
        self_conn = A(f, f);
        if useStandard
            % Find strong neighbors not in coarse set
            strongConn = part.strongConn(local);
            
            F_s = strongConn & ~P;
            
            if any(F_s)
                % Compute elimination coefficient and reduce old
                modifier = sum(N_conn(F_s)./D(N(F_s)));

                N_conn = N_conn./modifier;
                P_conn = P_conn./modifier;
                self_conn = self_conn./modifier;
                
                N_conn = N_conn(~F_s);
                
                % Add indirect neighbors via strong connections to vectors
                F_add = N(F_s);
                ignore = f;
                for j = 1:numel(F_add)
                    [xN, xP, xN_conn, xP_conn, xcoarseIndex] = getConn(F_add(j), ignore); %#ok
                    
                    N_conn      = [N_conn; xN_conn];%#ok
                    P_conn      = [P_conn; xP_conn];%#ok
                    coarseIndex = [coarseIndex; xcoarseIndex];%#ok
                    
                    ignore = [ignore; F_add(j)]; %#ok
                end
                
            end
        end        
        a = sum(N_conn)/sum(P_conn);
        vals = a*P_conn./self_conn;

        i_ind = [i_ind; repmat(f, numel(vals), 1)]; %#ok
        j_ind = [j_ind; coarseIndex]; %#ok
        v = [v; vals]; %#ok
    end

    prolong = sparse(i_ind, j_ind, v, size(A, 1), numel(coarse));
    prolong = bsxfun(@rdivide, prolong, sum(prolong, 2));
    
    restrict = prolong';
end

function [N, P, N_conn, P_conn, coarseIndex, local] = getConnections(f, ignore, rowi, coli, conn, coarse, part)
    if any(isnan(ignore))
        local = find(rowi == f);
        fun = @(x) any(x == f);
    else
        fun = @(x) any(x == f) && ~any(ismember(x, ignore));
        local = find(rowi == f & ~ismember(coli, ignore));
    end
    N = coli(local);

    % Neighboring connections - local mapping
    N_conn = conn(local);

    % Interpolatory set is all coarse blocks strongly connected to f
    coarseSubs = find(ismember(coarse, N));
    coarseIndex = find(cellfun(fun, part.strongCoarseConnected(coarseSubs)));
    coarseIndex = coarseSubs(coarseIndex); %#ok will be returned anyway
    
    P = ismember(N, coarse(coarseIndex));

    % Interpolatory connections - local mapping
    P_conn = N_conn(P);
end