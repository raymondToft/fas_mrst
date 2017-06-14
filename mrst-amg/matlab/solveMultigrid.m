function [x, res, meta, partitions] = solveMultigrid(G, A, q, varargin)
    opt = struct('iterations',  1000, ...
                 'tolerance',   1e-6, ...
                 'norm',        2, ...
                 'smoother',    'gauss-seidel', ...
                 'coarsening',  'uniform', ...
                 'interpolator','unitary', ...
                 'coarsefactor', 2, ...
                 'presmooth',    1, ...
                 'postsmooth',   1, ...
                 'levels',       G.griddim, ...
                 'cfstrength',   .25, ...
                 'fullweight',   true, ...
                 'maxCoarseSize',   5, ...
                 'maxCoarseLevels', inf, ...
                 'cycleindex',   1, ...
                 'cycletype',    [], ...
                 'jacobiBasis',  0, ...
                 'verbose',      mrstVerbose, ...
                 'partitions',   [], ...
                 'fullMG',       false, ...
                 'doSolve',      true, ...
                 'enforceVariational',  true, ...
                 'x0',           zeros(size(q)), ...
                 'omega',        2/3 ...
                 );
    opt = merge_options(opt, varargin{:});
    
    if opt.presmooth == 0 && opt.postsmooth == 0
        warning('Without smoothing, this will not converge.')
    end
    setup_timer = tic();
    
    if ~isempty(opt.cycletype)
        switch lower(opt.cycletype)
            case 'v'
                opt.cycleindex = 1;
            case 'w'
                opt.cycleindex = 2;
            case 'f'
                opt.cycleindex = [2 1];
            otherwise
                error('Unknown cycle!')
        end
    end
    switch opt.interpolator
        case 'mstpfa'
            warning('off','MATLAB:nearlySingularMatrix')
            warning('off','MATLAB:singularMatrix')
            getOperators = @(A, q, p) mstpfaInterpolator(A, q, p, opt.enforceVariational);
        case {'unitary', 'constant'}
            getOperators = @(A, q, p) unitInterpolator(A, q, p);
        case 'linear'
            getOperators = @(A, q, p) tensorLinearInterpolator(A, q, p, opt.fullweight);
        case 'standard'
            getOperators = @(A, q, p) multigridInterpolator(A, q, p, true);
        case 'direct'
            getOperators = @(A, q, p) multigridInterpolator(A, q, p, false);
        otherwise
            error('Unknown operator');
    end
    
    partitions = opt.partitions;
    if isempty(partitions)
        switch lower(opt.coarsening)
            case {'uniform', 'metis', 'pairwise'}
                partitions = makePartitions(G, A, q, getOperators, ...
                                            opt.coarsefactor, ...
                                            strcmpi(opt.interpolator, 'mstpfa'), ...
                                            opt);
            case 'standard'
                partitions = makeCoarseningFC(A, q, getOperators, opt);
                
            otherwise
                error()
        end
    end
    
    switch opt.smoother
        case 'jacobi'
            smoother = @(A, q, x, n) smootherJacobi(A, q, x, n, opt.omega);
        case {'gauss-seidel', 'gs'}
            smoother = @(A, q, x, n) smootherGS(A, q, x, n, true);
        case {'gauss-seidel-symm', 'gs-symm'}
            % Symmetric gauss-seidel - do a upper/lower followed by a
            % lower/upper pass.
            l1 = @(A, q, x, n) smootherGS(A, q, x, n, true);
            smoother = @(A, q, x, n) smootherGS(A, q, l1(A, q, x, n), n, false);
        otherwise
            error('Unknown smoother');
    end
    x   = opt.x0;
    res = zeros(opt.iterations+1, 1);
    res(1) = norm(A*x - q, opt.norm)/norm(q, opt.norm);
    
    meta = cell(opt.iterations, 1);
    
    if ~opt.doSolve
        return
    end
    
    solvemg = @(A, q, x, partitions, level) multigridCycle(A, q, x, ...
                                opt.presmooth, opt.postsmooth,...
                                opt.cycleindex, partitions, ...
                                smoother, getOperators, makeMeta(), level);
    if opt.fullMG
        levels = coarsenSystems(A, q, partitions, getOperators);
        defect = q - A*x;
    end
    dispif(opt.verbose, 'Spent %g seconds in setup phase!\n', toc(setup_timer));
    
    if opt.doSolve
        timer_solve = tic();
        for i = 1:opt.iterations
            if opt.fullMG
                % Full multigrid
                meta_tmp = makeMeta();
                % Project to the coarsest level
                x_local = projectDown(zeros(size(q)), partitions, getOperators);
                [defect, defects] = projectDown(defect, partitions, getOperators); %#ok
                for j = numel(partitions):-1:0
                    % Work ourselves up to the finest level, doing a multigrid cycle at
                    % each level
                    A_c = levels{j+1}.A;
                    if j == numel(partitions)
                        % Lowest level, solve exactly
                        x_local = A_c\defects{j+1};
                        meta_tmp.level  = [meta_tmp.level;  j];
                        meta_tmp.defect = [meta_tmp.defect; 0];
                    else
                        [x_local, m] = solvemg(A_c, defects{j+1}, x_local, partitions((j+1):end), j);
                        meta_tmp.level  = [meta_tmp.level;  m.level];
                        meta_tmp.defect = [meta_tmp.defect; m.defect];
                    end
                    if j ~= 0;
                        % We are not at the finest level, prepare data upwards
                        [R, I] = getOperators(A_c, defects{j+1}, partitions{j}); %#ok
                        % Interpolate to finer level
                        x_local = I*x_local;
                    end
                end
                meta{i} = meta_tmp;
                x = x + x_local;
                % Compute defect as right hand side for the next cycle
                defect = q - A*x;
            else
                % Repeated cycles until convergence
                [x, meta{i}] = solvemg(A, q, x, partitions, 0);
            end
            res(i+1) = norm(A*x - q, opt.norm)/norm(q, opt.norm);
            if res(i+1) < opt.tolerance
                dispif(opt.verbose, 'Converged in %d iterations\n', i);
                break
            end
            if res(i+1) >= res(i)
                warning('Stagnated in %d iterations\n', i);
                fprintf('Final residual: %1.4g - tolerance was %1.4g\n', res(i+1), opt.tolerance);
                break;
            end
        end
        if res(i+1) > opt.tolerance
            dispif(opt.verbose, 'Did not converge in %d iterations. Final residual: %1.4g\n', i, res(end));
        end

        dispif(opt.verbose, 'Spent %g seconds in solver phase!\n', toc(timer_solve));
    end
end

function partitions = addInterpolatorsSimple(A, rhs, partitions, getOperators)
    for i = 1:numel(partitions)
        G = partitions{i}.grid;
        n = prod(G.parent.cartDims);
        
        [R, I] = getOperators(A, rhs, partitions{i});
        partitions{i}.R = R;
        partitions{i}.I = I;
        
        A = R*A*I;
        rhs = R*rhs;
    end
end

function levels = coarsenSystems(A, q, partitions, getOperators)
    n = numel(partitions);
    levels = cell(n+1, 1);
    levels{1} = struct('A', A, 'q', q);
    for i = 1:n
        [R, I] = getOperators(A, q, partitions{i});
        A = R*A*I;
        q = R*q;
        levels{i+1} = struct('A', A, 'q', q);
    end
end

function [x, xs] = projectDown(x, partitions, getOperators)
    n = numel(partitions);
    xs = cell(numel(partitions)+1, 1);
    xs{1} = x;
    for i = 1:n
        % Assume caching
        [R, I] = getOperators([], [], partitions{i});
        x = R*x;
        xs{i+1} = x;
    end
end

function meta = makeMeta()
    meta = struct('level',  [], ...
                  'defect', []);
end
