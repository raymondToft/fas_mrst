 mrstModule add coarsegrid mex gridtools
%%
rng(0)
[G, A, q, rock] = getTestCase('dims', [2, 2].^3, ...
                            'perm', 'layers', ...
                            'layers', 3, ...
                             'gridaspect', 1);
figure(1); clf;
plotCellData(G, log10(rock.perm), 'edgec', 'w', 'edgea', .5)
axis tight off
%%
tic
p = coarsenColorizeAMG(A, 10, 0.2);
toc
%%
clf
plotDependence(G, p)