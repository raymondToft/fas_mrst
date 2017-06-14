function p = coarsenColorizeAMG(A, n, threshold)
    dependence = cell(n, 1);
    nc = size(A, 1);
    
    getIndicator = @(A) diag(A);
    indicator = getIndicator(A);
    C = zeros(nc, 1);
    
    iter = 1;
    categorized = false(nc, 1);
    
    while any(~categorized)
        
        [v, ind] = max(indicator);
        
        C(iter) = ind;
        dep = find(A(ind, :) ~= 0);
        
        indicator(dep) = 0;
        indicator(ind) = 0;
        
        categorized(dep) = 1;
        
        dependence{iter} = dep;
        iter = iter + 1;
    end
    C(C == 0) = [];
    
    p = struct('n', n, 'C', C);
    p.dependence = dependence;
end
