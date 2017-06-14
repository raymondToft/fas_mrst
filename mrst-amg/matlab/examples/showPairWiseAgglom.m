 
n = 42;
[G, A, q, rock] = getTestCase('dims', [n, n], 'perm', 'jump2d');
 
 
 close all
 
for i = 1:size(rock.perm, 2);
    figure
    plotCellData(G, rock.perm(:, i), 'EdgeColor', 'w')
    colorbar
    axis tight
end
%%
close all
Mp = max(rock.perm(:));
mp = min(rock.perm(:));

colors = parula(3);


h1 = plotGrid(G, rock.perm(:, 1) == Mp & rock.perm(:, 2) == mp, 'FaceColor', colors(:, 1));
h2 = plotGrid(G, rock.perm(:, 1) == mp & rock.perm(:, 2) == Mp, 'FaceColor', colors(:, 2));
h3 = plotGrid(G, rock.perm(:, 1) == Mp & rock.perm(:, 2) == Mp, 'FaceColor', colors(:, 3));
h4 = plotGrid(G, rock.perm(:, 1) == mp & rock.perm(:, 2) == mp, 'FaceColor', 'w');
legend([h1, h2, h3, h4], 'K_x=100, K_y=1', 'K_x=1, K_y=100', 'K_x=100, K_y=100', 'K_x=1, K_y=1');

axis tight
export_fig('../figures/grids/jump2d_perm', '-eps', '-pdf', '-transparent')
%%
solver = @(varargin) solveMultigrid(G, A, q,...
                                'coarsening', 'pairwise', ...
                                'interpolator', 'constant',...
                                'cfstrength',   .25, ...
                                'maxCoarseSize', 10, ...
                                'smoother', 'gs-symm', ...
                                'maxCoarseLevels', 8, ...
                                'levels',          1, ...
                                'coarsefactor',    cfactor, ...
                                'cycleindex',      2, ...
                                'verbose', true, varargin{:});
[x, res, meta, partitions] = solver();
%%
s = plotPartitions(G, partitions);
close all
for i = 1:2
    figure;
    h1 = plotGrid(G, rock.perm(:, 1) == Mp & rock.perm(:, 2) == mp, 'FaceColor', colors(:, 1), 'edgec', 'none');
    h2 = plotGrid(G, rock.perm(:, 1) == mp & rock.perm(:, 2) == Mp, 'FaceColor', colors(:, 2), 'edgec', 'none');
    h3 = plotGrid(G, rock.perm(:, 1) == Mp & rock.perm(:, 2) == Mp, 'FaceColor', colors(:, 3), 'edgec', 'none');
    outlineCoarseGrid(G, s.parts(:, i), 'k')
    axis tight
    export_fig(['../figures/grids/jump2d_grid_', num2str(i)], eps', '-pdf', '-transparent')

end