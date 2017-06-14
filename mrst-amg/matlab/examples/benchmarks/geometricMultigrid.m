name = 'cart';
name = 'l';
name = 'twist';
name = 'circle';
perm = 'homo';
perm = 'stripes';

[G, A, q, rock] = getTestCase('dims', [64 64], ...
                            'grid', name,...
                            'twistval', 0.01, ...
                            'perm', perm, ...
                            'layers', 3, ...
                            'stencil',  'tpfa', ...
                            'flow',   'bc');
                         
                         
figure(1); clf;
plotCellData(G, log10(rock.perm(:,1)), 'edgec', 'w', 'edgea', .5)
axis tight off
%%
figure(1); clf;
plotCellData(G, A\q)
axis tight 
fname = ['../figures/geometricmg/', perm, '_', name, '_reference'];
export_fig([fname, '.eps'], '-transparent')
export_fig([fname, '.pdf'], '-transparent')

%%
makeProblem = @(interpolator, coarsening, ...
                smoother, cycle) struct('interpolator', interpolator,...
                                                 'smoother',   smoother, ...
                                                 'coarsening', coarsening, ...
                                                 'cycle',  cycle);
                                                   
problems = {};

smoothers = {'jacobi', 'gauss-seidel'};
cycles = {'v', 'w', 'f'};
% p = {'linear', 'unitary'};
p = {'linear'};
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
    interpol = problems{i}.interpolator;
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
                                    'coarsefactor',    [8 2], ...
                                    'cycletype',       cycle, ...
                                    'tolerance',       1e-6, ...
                                    'iterations',      100, ...
                                    'verbose',         true);
    problems{i}.res = res;
    problems{i}.meta = meta;
    problems{i}.partitions = partitions;
    problems{i}.name = [interpol, ' ', smoother, ' (', cycle, '-cycle)'];
end


%%
close all
df = get(0, 'defaultfigureposition');
figure('position', df.*[1 1 1.25 .75]);


tmp = cellfun(@(x) x.res, problems, 'unif', false);
tmp = horzcat(tmp{:});
nonzero = sum(tmp, 2) > 0;
tmp = tmp(nonzero, :);


isLinear = cellfun(@(x) strcmpi(x.interpolator, 'linear'), problems);

semilogy(tmp(:, isLinear), '-', 'LineWidth', 2)
hold on
semilogy(tmp(:, ~isLinear), '--', 'LineWidth', 2)

ylabel('||Ax-b||_2/||b||_2')
xlabel('Iteration count')
% grid on
legend(cellfun(@(x) x.name, problems(isLinear), 'unif', false))
axis tight
xlim([0 100])
set(gcf,'renderer','painters')

fname = ['../figures/geometricmg/', perm, '_', name, '_convergence'];
% export_fig([fname, '.eps'], '-transparent')
% export_fig([fname, '.pdf'], '-transparent')

%%
close all
plotGrid(G, 'FaceColor', 'none')
outlineCoarseGrid(G, partitions{1}.partition, 'k')
axis tight off