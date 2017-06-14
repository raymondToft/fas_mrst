function [restrict, prolong] = unitInterpolator(A, q, p, varargin)
    cache = true;
    if numel(varargin)
        cache = varargin{1};
    end
    if isstruct(p)
        partition = p.partition;
    else
        partition = p;
    end
    
    if isfield(p, 'R') && isfield(p, 'I') && cache
        %caching for speedup
        restrict = p.R;
        prolong = p.I;
        return
    end
    
    n = numel(partition);
    m = max(partition);
    i = (1:n) .';
    j = partition(i);
    v = ones(n, 1);
    prolong = sparse(i, j, v, n, m);
    restrict = prolong .';
end