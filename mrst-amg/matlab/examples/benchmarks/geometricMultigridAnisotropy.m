% name = 'cart';
name = 'l';
% name = 'twist';
% name = 'circle';
% perm = 'homo';
% perm = 'stripes';
perm = 'anisotropy'

[G, A, q, rock] = getTestCase('dims', [64 64], ...
                            'grid', name,...
                            'twistval', 0.01, ...
                            'perm', perm, ...
                            'layers', 3, ...
                            'stencil',  'tpfa', ...
                            'flow',   'src');
                         
                         
figure(1); clf;
plotCellData(G, log10(rock.perm(:,1)), 'edgec', 'w', 'edgea', .5)
axis tight off
%%
close all
figure(1); clf;
plotCellData(G, A\q)
axis tight 
fname = ['../figures/geometricmg/', perm, '_', name, '_reference'];
export_fig([fname, '.eps'], '-transparent')
export_fig([fname, '.pdf'], '-transparent')
%%
makeProblem = @(coarsefactor, coarsening, ...
                smoother, cycle) struct('coarsefactor', coarsefactor,...
                                                 'smoother',   smoother, ...
                                                 'coarsening', coarsening, ...
                                                 'cycle',  cycle);
                                                   
problems = {};
smoothers = {'gauss-seidel'};
% smoothers = {'jacobi', 'gauss-seidel'};
cycles = {'v', 'f'};
p = {[8 2], [4 4], [2 8]};
for k = p
for i = 1:numel(smoothers)
    for j = 1:numel(cycles)
        problems = [problems; makeProblem(k,  'uniform', smoothers{i}, cycles{j})];
%         problems = [problems; makeProblem('unitary', 'uniform', smoothers{i}, cycles{j})];
    end
end
end

names = {};
for i = 1:numel(problems)
    fprintf('Solving %d of %d\n', i, numel(problems));
    coarse = problems{i}.coarsening;
    coarsefactor = problems{i}.coarsefactor;
    interpol = 'linear';
    smoother  = problems{i}.smoother;
    cycle = problems{i}.cycle;
    [x, res, meta, partitions] = solveMultigrid(G, A, q,...
                                    'coarsening',       coarse, ...
                                    'interpolator',     interpol,...
                                    'cfstrength',       .25, ...
                                    'maxCoarseSize',    8, ...
                                    'smoother',        smoother, ...
                                    'maxCoarseLevels', 5, ...
                                    'levels',          2, ...
                                    'coarsefactor',    coarsefactor, ...
                                    'cycletype',       cycle, ...
                                    'tolerance',       1e-6, ...
                                    'iterations',      200, ...
                                    'verbose',         true);
    problems{i}.res = res;
    problems{i}.meta = meta;
    problems{i}.partitions = partitions;
    problems{i}.name = [cycle, '-cycle Nx/Ny = ', num2str(coarsefactor(1)), '/', num2str(coarsefactor(2))];

%     problems{i}.name = [interpol, ' ', smoother, ' (', cycle, '-cycle), Nx/Ny = ', num2str(coarsefactor(1)), '/', num2str(coarsefactor(2))];
end


%%
close all
df = get(0, 'defaultfigureposition');
figure('position', df.*[1 1 1.25 .75]);


tmp = cellfun(@(x) x.res, problems, 'unif', false);
tmp = horzcat(tmp{:});
nonzero = sum(tmp, 2) > 0;
tmp = tmp(nonzero, :);



semilogy(tmp, '-', 'LineWidth', 2)

ylabel('||Ax-b||_2/||b||_2')
xlabel('Iteration count')
% grid on
legend(cellfun(@(x) x.name, problems, 'unif', false), 'location', 'southeast')
axis tight
xlim([0 200])
set(gcf,'renderer','painters')

fname = ['../figures/geometricmg/', perm, '_', name, '_convergence'];
export_fig([fname, '.eps'], '-transparent')
export_fig([fname, '.pdf'], '-transparent')

%%
close all

offset = (numel(smoothers)*numel(cycles));
colors = jet(5);
for i = 1:offset:(numel(problems)-offset+1)
    
    
    figure;
%     plotGrid(G, 'FaceColor', 'none')
    s = plotPartitions(G, problems{i}.partitions, false);
    close;
    figure;
    plotCellData(G, rock.perm(:,1), 'EdgeColor', 'k', 'EdgeAlpha', .2, 'LineSmoothing', 'on')
    colormap winter
    caxis([.75*min(rock.perm(:,1)), 1.25*max(rock.perm(:,1))])
    colormap jet
    outlineCoarseGrid(G, s.Partition_1, 'w')

%     outlineCoarseGrid(G, problems{i}.partitions{1}.partition, 'w')
    axis tight off
    
    coarsefactor = problems{i}.coarsefactor;
    fname = ['../figures/geometricmg/', perm, '_', name, num2str(coarsefactor(1)), '_', num2str(coarsefactor(2)),'_grid'];
    export_fig([fname, '.eps'], '-transparent')
    export_fig([fname, '.pdf'], '-transparent')
end
%%
close all
[G, A, q, rock] = getTestCase('dims', [8 8], ...
                            'grid', 'cart',...
                            'twistval', 0.01, ...
                            'perm', 'anisotropy', ...
                            'layers', 3, ...
                            'stencil',  'tpfa', ...
                            'flow',   'src');
plotCellVectorData(G, [5*ones(G.cells.num,1), zeros(G.cells.num, 1)], 'scale', .5, 'linewidth', 2, 'colormap', @(N) repmat([1 0 0], N, 1))
plotCellVectorData(G, [zeros(G.cells.num,1), ones(G.cells.num, 1)], 'scale', .1, 'linewidth', 2, 'colormap', @(N) repmat([0 0 1], N, 1))
plotGrid(G, 'facec', 'none')
axis tight off

%%
close all
T = computeTrans(G, rock);

facetrans = zeros(G.faces.num, 1);
facetrans(G.cells.faces(:, 1)) = T;

vtrans = bsxfun(@(x,y) x.*y, G.faces.normals, facetrans);
vtrans = [vtrans; -vtrans];
X = G.faces.centroids(:, 1);
Y = G.faces.centroids(:, 2);
ptypes = {[2 2], [1 4], [4 1]};
for i = 1:numel(ptypes)
    figure;
    h = quiver([X; X], [Y; Y], vtrans(:, 1), vtrans(:, 2), .5, 'LineWidth', 2, 'MaxHeadSize', 1, 'color', 'r');
    
    p = partitionUI(G, ptypes{i});
    outlineCoarseGrid(G, p, 'b')
    plotGrid(G, 'facec', 'none')
    axis tight off
    
    
    
    fname = ['../figures/geometricmg/', name, '_quiver', num2str(i)];
    export_fig([fname, '.eps'], '-transparent')
    export_fig([fname, '.pdf'], '-transparent')
end
