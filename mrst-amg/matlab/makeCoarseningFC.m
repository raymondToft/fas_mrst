function partitions = makeCoarseningFC(A, q, getOperators, opt)

    partitions = {};
    coarseNo = inf;
    iter = 1;
    while coarseNo > opt.maxCoarseSize  && iter <= opt.maxCoarseLevels
        part = standardCoarsening(A, opt.cfstrength);
        partitions = [partitions; part]; %#ok
        
        [R, I] = getOperators(A, q, partitions{end});
        if opt.jacobiBasis
            I = jacobiUpdateBasis(A, I, opt.jacobiBasis, 2/3, opt.cfstrength);
            I = bsxfun(@rdivide, I, sum(I, 2));
            R = I';
        end
        partitions{end}.R = R;
        partitions{end}.I = I;

        A = R*A*I; q = R*q;
        
        partitions{end}.A = A;
        partitions{end}.q = q;
        
        coarseNo = sum(part.coarse);
        
        dispif(opt.verbose, 'Partitioned level %d of maximum %d\n', iter, opt.maxCoarseLevels);
        iter = iter + 1;
    end
end

