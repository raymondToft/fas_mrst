 [G, A, q, rock] = getTestCase('dims', [128, 1], ...
                              'grid', 'cart',...
                              'perm', 'homo', ...
                              'stencil',  'tpfa', ...
                              'flow',   'none', ...
                              'griddims',  [1 1], ...
                              'ordering', 'default');
                          
                          
 [~, A_coarse, q_coarse, ~] = getTestCase('dims', [32, 1], ...
                              'grid', 'cart',...
                              'perm', 'homo', ...
                              'stencil',  'tpfa', ...
                              'flow',   'none', ...
                              'griddims',  [1 1], ...
                              'ordering', 'default');
gs = @(x, n) smootherGS(A, q, x, n);

X = G.cells.centroids(:, 1);
getData = @(K) sin(X*pi*K);
%%
close all
x0 = getData(60) + getData(8);
df = get(0, 'defaultfigureposition');
figure('position', df.*[1 1 1.8 .75]);


subplot(2,2,[1 3])

x = x0;
x = smootherGS(A, q, x, 10);
hold on
% plot(A*x - q, 'r')
plot(x0, '--', 'LineWidth', 2)
plot(x, 'r--', 'LineWidth', 2)
plot(0,0, 'g--', 'LineWidth', 2')
grid on
axis tight
title('First step, fine grid')
legend('Initial', 'First iterate', 'Second iterate', 'location', 'south',...
                                'orientation', 'horizontal')

subplot(2, 2, 2)

hold on

x_coarse0 = x(1:4:(numel(q)-3));
plot(x_coarse0, 'ro--', 'LineWidth', 2)
x_coarse = smootherGS(A_coarse, q_coarse, x_coarse0, 10);

plot(x_coarse, 'go--', 'LineWidth', 2)
grid on

axis tight
ylim([min(x0), max(x0)])
title('Second step, coarse grid')


subplot(2, 2, 4)
hold on
plot(x, 'r--', 'LineWidth', 2)
xnew = smootherGS(A, q, x, 10);
plot(xnew, 'g--', 'LineWidth', 2)
grid on
axis tight
ylim([min(x0), max(x0)])
title('Second step, fine grid')
% export_fig(['../figures/smoothers/aliasing.eps'], '-transparent');