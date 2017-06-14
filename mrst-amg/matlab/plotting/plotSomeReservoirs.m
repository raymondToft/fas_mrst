%%
load /data/gullfaks

%%
close all;
df = get(0, 'defaultfigureposition');
figure('position', df.*[1 1 1.75 1.25]);



h = plotCellData(G, log10(rock.perm(:,1)), 'EdgeColor', [1 1 1]*.2, 'LineWidth', .2);
view(-135, 32)
c = colorbar();
axis tight
box on

% tmp = get(c, 'TickLabels');
% set(c, 'TickLabels', cellfun(@(x) ['10^{', x, '}'], tmp, 'UniformOutput', false))
% parula
colormap jet
fname = ['../figures/grids/gullfaks'];
export_fig([fname, '.png'], '-transparent', '-m2')
% print(gcf, '-dpng', 'fname')

%%
load /data/norne

%%
close all;
df = get(0, 'defaultfigureposition');
figure('position', df.*[1 1 1.75 1.25]);



h = plotCellData(G, log10(rock.perm(:,1)), 'EdgeColor', [1 1 1]*.2, 'LineWidth', .2);
view(77, 40)
c = colorbar();
axis tight
box on

% tmp = get(c, 'TickLabels');
% set(c, 'TickLabels', cellfun(@(x) ['10^{', x, '}'], tmp, 'UniformOutput', false))
% parula
% colormap parula
colormap jet
fname = ['../figures/grids/norne'];

export_fig([fname, '.png'], '-transparent', '-m2')
% print(gcf, '-dpng', 'fname')
%%
%%
mrstModule add spe10 mex
[G, w, rock] = SPE10_setup();
close all;
df = get(0, 'defaultfigureposition');
figure('position', df.*[1 1 1.75 1.25]);



h = plotCellData(G, log10(rock.perm(:,1)), 'EdgeColor', [1 1 1]*.2, 'LineWidth', .2);
view(77, 40)
c = colorbar();
axis tight
box on

% tmp = get(c, 'TickLabels');
% set(c, 'TickLabels', cellfun(@(x) ['10^{', x, '}'], tmp, 'UniformOutput', false))
% parula
% colormap parula
colormap jet
fname = ['../figures/grids/spe10'];

export_fig([fname, '.png'], '-transparent', '-m2')
% print(gcf, '-dpng', 'fname')