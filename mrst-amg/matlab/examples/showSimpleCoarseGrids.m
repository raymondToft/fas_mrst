 [G, A, q, rock] = getTestCase('dims', [64, 64], ...
                              'grid', 'cart',...
                              'perm', 'homo', ...
                              'stencil',  'tpfa', ...
                              'flow',   'none', ...
                              'griddims',  [1 1], ...
                              'ordering', 'default');
%%
[x, res, meta, partitions] = solveMultigrid(G, A, q, 'doSolve', false, 'maxCoarseSize', 0);
%%
close all
s = plotPartitions(G, partitions);
close all
plotToolbar(G, struct('s', s, 'rock', rock))
axis tight
%%

close all

np = numel(partitions);
colors = jet(np+1);

for pi = 1:(np+1)
    
    
    figure
    for i = 1:pi;
        if i == 0
             plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [1 1 1].*.6)
        elseif i <= np
            outlineCoarseGrid(G, s.parts(:, i), colors(i, :))
        else
            outlineCoarseGrid(G, ones(G.cells.num, 1), colors(i, :))
        end
        axis off
    end
%     export_fig(['../figures/grids/simple_coarsening', num2str(pi), '.eps'], '-transparent');
end

%% Show interpolation stuff
n = 30;
 [G, A, q, rock] = getTestCase('dims', [n 1], 'griddims',  [1 1]);
M = 10;
 
part = repmat((1:M), n/M, 1);
part = part(:);
part = processPartition(G, part);
part = compressPartition(part);

plotCellData(G, part);

CG = generateCoarseGrid(G, part);
CG = coarsenGeometry(CG);
partition = struct('grid', CG, 'centers', mapCenters(CG, G), 'partition', part);

[restrict, prolong] = tensorLinearInterpolator(A, q, partition, true);

[~, Ic] = unitInterpolator(A, q, part);
close all

plotToolbar(G, restrict')
%%
if 1
    fun = @(x) sin(x*pi);
    name = 'sine';
    lpos = 'northeast';
else
    fun = @(x) double(x > 0.5);
    name = 'delta';
    lpos = 'northwest';
end
close all
df = get(0, 'defaultfigureposition');
figure('position', df.*[1 1 1.8 .75]);

x_f = G.cells.centroids(:, 1);
x_c = G.cells.centroids(partition.centers, 1);



f_c = fun(x_c);
f_f = fun(x_f);
f_i_linear = prolong*f_c;
f_i_const = Ic*f_c;


hold on
area(x_f, f_f, 'FaceColor', [.3 .3 1])

plot(x_f, f_i_linear, 'ok', 'MarkerFaceColor', 'y')
plot(x_f, f_i_const, 'ok', 'MarkerFaceColor', 'b')

plot(x_c, f_c, 'ok', 'MarkerFaceColor', 'r')

grid on
legend({'f(x) on fine', 'f(x) linear interpolation', 'f(x) constant interpolation', 'f(x) on coarse'}, 'location', lpos)
export_fig(['../figures/grids/1d_interp_', name ,'.eps'], '-transparent')
%%
close all
hold on
for i = 1:size(prolong, 2)
    plot(prolong(:, i))
end
%%
% name = 'cart';
% name = 'l';
% name = 'twist';
name = 'circle';
dims = [32, 32];
[G, A, q, rock] = getTestCase('dims', dims, ...
                              'grid', name,...
                              'perm', 'homo', ...
                              'stencil',  'tpfa', ...
                              'flow',   'bc', ...
                              'griddims',  [1 1], ...
                              'ordering', 'default');

G_embed = getTestCase('dims', dims, ...
                              'grid', 'cart',...
                              'perm', 'homo', ...
                              'stencil',  'tpfa', ...
                              'flow',   'bc', ...
                              'griddims',  [1 1], ...
                              'ordering', 'default');
%%
rehash
[x, res, meta, partitions] = solveMultigrid(G, A, q, 'doSolve', true, ...
                                'maxCoarseSize', 5, ...
                                'coarsefactor', 2, ...
                                'interpolator', 'linear',...
                                'verbose', true);

%%
s = plotPartitions(G, partitions, false);
close all
plotGrid(G, 'FaceColor', 'none')
outlineCoarseGrid(G, s.parts(:, 1),'r')
outlineCoarseGrid(G, s.parts(:, 2),'g')
outlineCoarseGrid(G, s.parts(:, 3),'b')
axis tight off
% export_fig(['../figures/grids/', name, '_partitioned.png'], '-transparent')
%%

close all
plotGrid(G)
% plotGrid(G_embed, 'FaceColor', 'None')
axis off
% export_fig(['../figures/grids/', name, '_embedded.png'], '-transparent')