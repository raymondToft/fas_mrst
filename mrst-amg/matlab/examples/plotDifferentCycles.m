 [G, A, q, rock] = getTestCase('dims', [16 16], ...
                            'perm', 'homo', ...
                            'layers', 3, ...
                            'flow', 'wells', ...
                             'gridaspect', [1 1]);
 solve = @(cindex, useFullMG, varargin) solveMultigrid(G, A, q, ...
                                            'tolerance',         1e-6,...
                                            'cycleindex',        cindex,...
                                            'presmooth',         1, ...
                                            'fullMG',            useFullMG, ...
                                            'postsmooth',        1, ...
                                            'coarsefactor',      2, ...
                                            'maxCoarseLevels', 3,...
                                            varargin{:});
saveplot = @(name) export_fig(['../figures/cycles/', name], '-transparent', '-eps', '-pdf');
df = get(0, 'defaultfigureposition');
fig = @() figure('position', df.*[1 1 1.25 1]);
%%

close all
fig()
[x, res, meta] = solve(1, false);
plotLevels(meta{1}.level)
saveplot('v')

fig()
[x, res, meta] = solve(2, false);
plotLevels(meta{1}.level)
saveplot('w')

fig()
[x, res, meta] = solve(3, false);
plotLevels(meta{1}.level)
saveplot('wv')

fig()
[x, res, meta] = solve([2 1], false);
plotLevels(meta{1}.level)
area([2 9], [2 2], 'FaceColor', c, 'Edgecolor', 'r', 'Annotation', 'off', 'LineStyle', '--', 'LineWidth', 2)
area([10 14], [2 2], 'FaceColor', c, 'Edgecolor', 'g', 'Annotation', 'off', 'LineStyle', '--', 'LineWidth', 2)

legend('Smoother', 'Exact solver', 'Embbedded w-cycle', 'Embbedded v-cycle')
saveplot('f')


%%
c = 'none';
ec = 'r';

close all
% c = [1 1 1]*.3;
for i = 1:2

    fig()
    [x, res, meta] = solve(i, true);
    hold on
    plotLevels(meta{1}.level)
    area([2 4], [1 1], 'FaceColor', c, 'Edgecolor', ec, 'Annotation', 'off', 'LineStyle', '--', 'LineWidth', 2)
    if i == 1
        a = [5 9]; b = [10 16];
    else
        a = [5 12]; b = [13 30];
    end
    area(a, [2 2], 'FaceColor', c, 'Edgecolor', ec, 'Annotation', 'off', 'LineStyle', '--', 'LineWidth', 2)
    area(b, [3 3], 'FaceColor', c, 'Edgecolor', ec, 'Annotation', 'off', 'LineStyle', '--', 'LineWidth', 2)

    legend('Smoother', 'Exact solver', 'Subcycle', 'Location', 'northwest')
    saveplot(['fullmg_', num2str(i)])

end

%%
figure;
[x, res, meta] = solve(2, true);
plotLevels(meta{1}.level)
