mrstModule add gridtools mrst-gui coarsegrid
%%
[G, A, q, rock] = getTestCase('dims', [60 60 1], ...
                            'grid', 'norne',...
                            'twistval', 0.01, ...
                            'perm', 'grid', ...
                            'layers', 3, ...
                            'stencil',  'tpfa', ...
                            'flow',   'bc', ...
                            'griddims',  [], ...
                            'ordering', 'default', ...
                            'gridaspect', [2 1 1/10]);
figure(1); clf;
plotToolbar(G, log10(rock.perm), 'edgec', 'w', 'edgea', .5)
axis tight off
figure(2); clf;
plotCellData(G, A\q, 'edgec', 'w', 'edgea', .5)

%%
maxCoarse = 200;
cfactor = 10;
if 1
    interpol = 'unitary';
    coarse = 'pairwise';
elseif 0
    interpol = 'mstpfa';
    coarse = 'metis';
else
    interpol = 'linear';
    coarse = 'uniform';
    maxCoarse = 10;
    cfactor = 2;
end
    maxCoarse = 10;
    cfactor = 2;
solver = @(varargin) solveMultigrid(G, A, q,...
                                'coarsening', coarse, ...
                                'interpolator', interpol,...
                                'cfstrength',   .001, ...
                                'maxCoarseSize', maxCoarse, ...
                                'smoother', 'gs', ...
                                'fullmg', true, ...
                                'cycletype', 'f', ...
                                'maxCoarseLevels', 8, ...
                                'maxCoarseSize', 400, ...
...%                                 'partitions',      partitions, ...
                                'levels',          2, ...
                                'coarsefactor',    cfactor, ...
                                'cycleindex',      2, ...
                                'verbose', true, varargin{:});
[x, res, meta, partitions] = solver();
%%
close all
s = plotPartitions(G, partitions);
close all
plotToolbar(G, struct('s', s, 'rock', rock))
axis tight
%%
i = 2;
close all
plotToolbar(G, log10(rock.perm))
outlineCoarseGrid(G, s.(['Partition_', num2str(i)]), 'k')

axis off equal

%%
for i = 1:numel(partitions)
    tmp = sum(partitions{i}.I, 2);
    full([min(tmp), max(tmp)])
end

%%
ind = 2;
A = partitions{ind}.A;

tic
part = partitionAgglomerationPairwise(A, 0.001, false);
toc
%%
clc
for i = 1:size(A, 2)
    a = A(i, :);
    a = a(a < 0);
    fprintf('%1.1g\t', full(a));
    fprintf('\n')
    
end
%%
tic
part = partitionAgglomerationLevels(A, 0.25, true, 3);
toc
%%
part = partitionAgglomerationLevels(A, .25, true, 5);