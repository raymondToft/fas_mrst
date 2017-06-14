mrstModule add coarsegrid mex gridtools mrst-gui
%%
% ga = [2 1] for spe10
[G, A, q, rock] = getTestCase('dims', [60 60], ...
                            'grid', 'cart',...
                            'twistval', 0.01, ...
                            'perm', 'tarbert', ...
                            'layers', 3, ...
                            'speoffset', 1, ...
                            'stencil',  'tpfa', ...
                            'flow',   'bc', ...
                             'gridaspect', [2 1]);
figure(1); clf;
plotCellData(G, log10(rock.perm(:, 1)), 'edgec', 'w', 'edgea', .5)
axis tight off
ref =  A\q;
figure(2); clf;
plotCellData(G, ref, 'edgec', 'w', 'edgea', .5)

%%
tic
part = standardCoarsening(A, .25);
toc
%%
close all
plotGrid(G, part.coarse, 'FaceAlpha', .5, 'facec', 'w')
plotCellData(G, log10(rock.perm))
%%
clf
plotGrid(G, part.coarse)

%%
L = 2;
P = 1;
p = coarseningByPath(A, G, part, L, P);

clf

plotCellData(G, log10(rock.perm))
plotGrid(G, p.coarse)

%%
tmp = p;
% tmp = part;
C = find(tmp.coarse);
for i = 1:numel(tmp.strongCoarseConnected)
    figure(1); clf;
    plotGrid(G, 'facec', 'none')
    plotGrid(G, tmp.strongCoarseConnected{i})
    plotGrid(G, C(i), 'FaceColor', 'r')
    pause
end
%%
tic
[restrict, prolong] = multigridInterpolator(A, q, part, true);
toc
%%
close all
plotToolbar(G, prolong, 'FaceAlpha', 1, 'EdgeColor', 'w', 'LineWidth', 2)
plotCellData(G, rock.poro./max(rock.poro))
axis tight

%%
% interpol = 'direct';
if 1
    interpol = 'standard';
    coarse = 'standard';
else
    interpol = 'linear';
    coarse = 'uniform';
end
[x, res, meta, partitions] = solveMultigrid(G, A, q,...
                                'coarsening', coarse, ...
                                'interpolator', interpol,...
                                'cfstrength',   .25, ...
                                'fullMG',       true, ...
                                'maxCoarseLevels', 3,...
                                'maxCoarseSize', 200, ...
                                'verbose', true);

%%
plotter = 'circle';
% plotter = 'simple';
close all
% plotCellData(G, numel(partitions)*(rock.poro./max(rock.poro)))
visualizeCFSplit(G, partitions, plotter)
axis off tight
%%
d = visualizeInterpolators(G, partitions);
close all

plotToolbar(G, d, 'FaceAlpha', 1)
axis tight off
figure;
plotCellData(G, log10(rock.perm), 'edgec', 'w', 'edgea', .5)
axis tight off
%%
close all
subplot(1,2,1)
plotToolbar(G, mod(d.level_4.partition, 13))
axis equal tight off

subplot(1,2,2)
plotCellData(G, log10(rock.perm(:,1)))
axis equal tight off
%%
close all
outlineCoarseGrid(G, d.level_5.partition, 'k')
plotCellData(G, log10(rock.perm(:,1)))
axis equal tight off

%%
for i = 1:numel(partitions)
    I = partitions{i}.I;
    d = sum(I, 2);
    [min(d), max(d)]
%     min(I(:))
end