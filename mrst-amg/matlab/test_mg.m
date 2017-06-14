mrstModule add coarsegrid mex gridtools mrst-gui
%%
[G, A, q, rock] = getTestCase('dims', [16 16], ...
                            'perm', 'homo', ...
                            'layers', 3, ...
                            'flow', 'wells', ...
                             'gridaspect', [1 1]);
figure(1); clf;
plotCellData(G, log10(rock.perm), 'edgec', 'w', 'edgea', .5)
axis tight off
ref =  A\q;
%%
p1 = partitionUI(G, ceil(G.cartDims./5));

CG = generateCoarseGrid(G, p1);
CG = coarsenGeometry(CG);
% Make basic interpolator
p = struct('partition', p1, 'grid', CG);
[R, I] = unitInterpolator(A, q, p);
subplot(2,2,1); spy(A);
subplot(2,2,2); spy(R*A*I);
subplot(2,2,3); plotCellData(G, ref);
axis tight off
coarse = R*A*I \ R*q;
subplot(2,2,4); plotCellData(G, coarse(p1));
axis tight off
%%
% interp = 'unitary';

if 1
    interp = 'linear';
    coarsening = 'uniform';
else
    interp = 'mstpfa';
%     coarsening = 'metis';
    coarsening = 'uniform';
end

partitions = [];

% cyc = 3;

% cyc = @(level) level;
smoother = 'gauss-seidel';
solve = @(varargin) solveMultigrid(G, A, q, 'tolerance',         1e-6,...
                                            'cycleindex',        1,...
                                            'presmooth',         1, ...
                                            'fullMG',            false, ...
                                            'cyclefun',          @(level) [], ...
                                            'postsmooth',        1, ...
                                            'coarsefactor',      2, ...
                                            'verbose',           true, ...
                                            'coarsening',        coarsening, ...
                                            'partitions',        partitions, ...
                                            'interpolator',      interp, ...
                                            varargin{:});

[x, res, meta, partitions] = solve( 'smoother',  smoother);


%%
close all
plotToolbar(G, struct('mg', x, 'mldivide', ref))
figure;
semilogy(res, 'linewidth', 2)
set(gcf, 'renderer', 'painters')

%%
close all
plotMultigridConvergence(meta)
%%

close all
plotPartitions(G, partitions)
%%
% [restrict, prolong] = tensorLinearInterpolator([], [], partitions{1}, true, false);
prolong = partitions{1}.I;
plotToolbar(G, prolong)

[min(sum(prolong, 2)), max(sum(prolong, 2))]

% outlineCoarseGrid(G, partitions{1}.partition)
colorbar
%%
index = 1;
coarsegrids = cellfun(@(x) x.grid, partitions, 'UniformOutput', false);

B = partitions{index}.I;
for i = 1:size(B, 2)
    figure(1); clf;
    data = B(:, i);
    plotCoarseData(coarsegrids, index-1, data)
    plotCoarseData(coarsegrids, index-1, zeros(size(B, 1), 1), 'FaceColor', 'none', 'EdgeColor', 'k', 'EdgeAlpha', .2)
    pause
end
%%
standardCoarsening(A, .25)

