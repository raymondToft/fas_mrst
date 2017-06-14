function amgGridTesting(interpolator, c, basisSmooth, solvetype)
if 0
    interpolator =  'standard';
    interpolator = 'direct';
    interpolator = 'constant';
    c = 4;
end
aspect = [1 1];
levels = 10;
mc = 10;
speoffset = 0;
cfs = .25;
switch c
    case 1
        dims = [32 32];
        perm = 'homo';
        grid = 'cart';
    case 2
        dims = [32 32];
        perm = 'anisotropy';
        grid = 'cart';
    case 3
        dims = [60 220];
        perm = 'tarbert';
        grid = 'cart';
        aspect = [1200/60, 2200/220]*ft;
    case 4
        dims = [60 220];
        perm = 'ness';
        grid = 'cart';
        aspect = [1200/60, 2200/220]*ft;
        speoffset = 10;
    case 5
        dims = [1 1];
        perm = 'grid';
        grid = 'norne';
        cfs = .001;
        mc = 500;
    case 6
        dims = [1 1];
        perm = 'grid';
        grid = 'gullfaks';
        %cfs = .001;
        mc = 500;
    case 7
        dims = [30, 30, 10];
        aspect = [10 10 1];
        perm = 'layers';
    otherwise
        warning('Bad case')
        return
end

[G, A, q, rock] = getTestCase('dims', dims, ...
                            'grid', grid,...
                            'twistval', 0.1, ...
                            'perm', perm, ...
                            'layers', 10, ...
                            'speoffset', speoffset, ...
                            'stencil',  'tpfa', ...
                            'flow',   'bc', ...
                             'gridaspect', aspect);
if basisSmooth
    tmp = ['_j', num2str(basisSmooth)];
else
    tmp = '';
end


if 0
    ratio = @(v) max(v(:))/min(v(:));
    str = [perm, ' & ', num2str(condest(A),3), ' & ', num2str(G.cells.num), ' & ' ...
           num2str(max(rock.perm(:)),3), ' & ',...
           num2str(min(rock.perm(:)),3), ' & ',...
           num2str(max(G.faces.areas),3), ' & ',...
           num2str(min(G.faces.areas),3),  '\\\\\n'];
%            num2str(ratio(rock.perm),3), ' & ', num2str(ratio(G.faces.area),3), '\\\\\n'];
       
    fprintf(str)
    return
end
% close all
% plotCellData(G, log10(rock.perm(:,1)));
% return
%% Initial solve to get partitions
isAgglom = strcmpi(interpolator, 'constant');

if isAgglom
    coarse = 'pairwise';
else
    coarse = 'standard';
end


isPreconditioned = false;
solver = @(varargin) solveMultigrid(G, A, q,...
                                'cfstrength',       cfs, ...
                                'maxCoarseLevels',  levels,...
                                'coarsening',       coarse, ...
                                'iterations',       100, ...
                                'maxCoarseSize',    mc, ...
                                'jacobiBasis',      basisSmooth, ...
                                'interpolator',     interpolator, ...
                                'verbose',          true,...
                                varargin{:});


cachename = fullfile(pwd, [grid, '_', perm, '_', coarse, '_', interpolator, '_', num2str(basisSmooth), '_partitions.mat']);

if exist(cachename, 'file')
    load(cachename, 'partitions');
else
    [x, res, meta, partitions] = solver('doSolve', false);
    save(cachename, 'partitions', '-v7.3');
end
    
solver = @(varargin) solver('partitions', partitions, varargin{:});

if ~isempty(solvetype)
    isPreconditioned = true;
    getPreconditioner = @(varargin) makeMultigridPreconditioner(G, A, q, partitions,...
                                    'cfstrength',       cfs, ...
                                    'maxCoarseLevels',  10,...
                                    'iterations',       1, ...
                                    'coarsening',       coarse, ...
                                    'maxCoarseSize',    1, ...
                                    'jacobiBasis',      basisSmooth, ...
                                    'interpolator',     interpolator, ...
                                    'verbose',          false, varargin{:});
    switch lower(solvetype)
        case 'gmres'
%             tol = 1e-7;
%             if strcmpi(grid, 'gullfaks')
                tol = 1e-9;
%             end
            solver = @(preconditioner) gmres(A, q, [], tol, 100, preconditioner);
        otherwise
            error('unknown solver')
    end
    tmp = [tmp, '_', lower(solvetype)];
end


savefig = @(name, varargin) export_fig(['../figures/amg/', grid, '_', perm, '_', interpolator, tmp, '_', name], ...
                                                        '-transparent', varargin{:});

%%
packProblem = @(cycle, fullmg, smoother) struct('fullMG', fullmg, 'cycletype', cycle, 'smoother', smoother);

problems = {};

cycles = {'v', 'f'};
smoothers = {'jacobi', 'gauss-seidel'};

for fullmg = [false, true]
    for k = smoothers
        for i = cycles            
            problems = [problems; packProblem(i{1}, fullmg, k{1})];
        end
    end
end


names = {};

for i = 1:numel(problems)
    fprintf('Solving problem %d of %d\n', i, numel(problems));

    p = problems{i};
    argin = {};
    for j = fieldnames(p)'
        argin = [argin; j; p.(j{1})];
    end
    if isPreconditioned
        precond = getPreconditioner(argin{:});
        res = zeros(1000, 1);
        
        [x, flag, relres, iter, ress] = solver(precond);
        res(1:numel(ress)) = ress;%./max(norm(q, 2), 1);
    else
        [x, res] = solver(argin{:});
    end
    
    problems{i}.res = res;

    name = '';
    if p.fullMG
        name = [name, ' ', 'FMG'];
    end
    name = [name, ' ', p.cycletype, '-cycle ', p.smoother];
    names = [names; name];
end
%%

close all
residuals = cellfun(@(x) x.res, problems, 'unif', false);
residuals = horzcat(residuals{:});
mid = numel(problems)/2;

semilogy(residuals(:, 1:mid), '--', 'LineWidth', 2)
hold on
semilogy(residuals(:, (mid+1) : end), '-', 'LineWidth', 2)
set(gcf, 'Renderer', 'painters')
legend(names)
ylim([1e-6 1])

oldlim = xlim();
xlim([1 min(oldlim(2), 100)])

savefig('convergence', '-eps', '-pdf')
% axis tight

if isPreconditioned
    return
end
%%
if ~isAgglom
    plotter = 'circle';
    % plotter = 'simple';
    close all


    if G.cells.num > 4096
        plotCellData(G, log10(rock.perm(:,1)))
        visualizeCFSplit(G, partitions, plotter, 3:numel(partitions));
    else
        plotGrid(G, 'FaceColor', 'none')
        visualizeCFSplit(G, partitions, plotter);
        h = colorbar();
        caxis([1, numel(partitions)])
        set(h, 'YTick', 1:6)
    end

    levels = visualizeCFSplit(G, partitions, plotter);
else
    close all
    s = plotPartitions(G, partitions);
    np = numel(partitions);
    colors = lines(np);
    for i = 1:2:np
        figure;
        if G.griddim == 2
            plotToolbar(G, log10(rock.perm(:,1)))
            h = outlineCoarseGrid(G, s.parts(:, i), 'b');
            if i < np
                h = outlineCoarseGrid(G, s.parts(:, i+1), 'r');
            end
            axis off equal
            set(h, 'LineWidth', 1)
        end
%         if G.griddim == 2
%             set(h, 'LineWidth', i)
%         end
        if G.griddim < 3
            savefig(['coarselevels_', num2str(i)], '-eps', '-pdf')
        end
    end
    
end

%%
axis equal off tight
if G.griddim < 3
savefig('coarselevels', '-eps', '-pdf')
end
%%
figure;
plotCellData(G, log10(rock.perm(:,1)))
axis equal off tight
savefig('permeability', '-eps', '-pdf')
%%
if isAgglom
    iter = 1:numel(partitions);
else
    iter = 1:2:numel(partitions);
end
if G.griddim == 3
    close all
    for i = iter
        figure
        if isAgglom
            plotCellData(G, mod(s.parts(:, i), 13));
        else
            plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', .05);
            plotCellData(G, ones(G.cells.num, 1)*i, levels == i, 'FaceAlpha', 1);
            plotCellData(G, (ones(G.cells.num, 1)*i+1), levels == i+1, 'FaceAlpha', 1);
            caxis([1 numel(partitions)]);
        end
        axis tight off
    %     camlight right
        if strcmpi(grid, 'norne')
            view(76, 43)
        elseif strcmpi(grid, 'gullfaks')
            view(-135, 52)
        end
        savefig(['coarselevel3d_', num2str(i)], '-png', '-m2')
    end
end
%%
close all
