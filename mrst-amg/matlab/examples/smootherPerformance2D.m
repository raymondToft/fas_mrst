[G, A, q, rock] = getTestCase('dims', [100 100], ...
                              'grid', 'cart',...
                              'perm', 'homo', ...
                              'stencil',  'tpfa', ...
                              'flow',   'none', ...
                              'griddims',  [1 1], ...
                              'ordering', 'default');
                        
                        
gs = @(x, n) smootherGS(A, q, x, n);

saveplot = @(name) export_fig(['../figures/smoothers/2d_', name, '.png'], '-transparent');

%%
close all
rng(0)


X = G.cells.centroids(:, 1);
Y = G.cells.centroids(:, 2);

plotter = @(x) surf(reshape(X, G.cartDims),...
                    reshape(Y, G.cartDims),...
                    reshape(x, G.cartDims));

K1 = 1;
K2 = 10;

x0 = sin(X*pi*K1) + sin(Y*pi*K2) + rand(G.cells.num, 1);

plotter(x0);
title('Initial data')
shading flat
zlim([-0.5 3])
saveplot('initial_surf')
x = x0;
for i = 1:1000
    x = smootherGS(A, q, x, 1);
    
    if ismember(i, [1 10 20 50 500 1000])
        figure
        plotter(x);
        shading flat
        zlim([-0.5 3])
        caxis([min(x0), max(x0)])
        title(['After ' num2str(i) ' Gauss-Seidel passes']);
        saveplot([num2str(i), '_surf'])
    end
end