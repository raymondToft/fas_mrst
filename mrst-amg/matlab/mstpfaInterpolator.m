function [restrict, prolong] = mstpfaInterpolator(A, q, p, useFE, varargin)
    cache = true;
    if numel(varargin)
        cache = varargin{1};
    end
    
    if isfield(p, 'R') && isfield(p, 'I') && cache
        %caching for speedup
        restrict = p.R;
        prolong = p.I;
        return
    end
    
    A = A - diag(sum(A, 2));
    CG = p.grid;
    
    basis = createFaceBasis(CG, A);
    [prolong, restrict] = assembleCoarseOperatorsPartition(CG, basis);
    if useFE
        restrict = prolong';
    end
end