 mrstModule add coarsegrid mex gridtools
%%

rng(0)
[G, A, q, rock] = getTestCase('dims', [3, 3].^3, ...
                              'perm', 'homo', ...
                              'flow', 'none', ...
                              'layers', 3, ...
                              'gridaspect', 1);
figure(1); clf;
plotCellData(G, log10(rock.perm), 'edgec', 'w', 'edgea', .5)
axis tight off
%%
solvers = {};
pack = @(fun, name) struct('f', fun, 'name', name);


jacobi = @(x) smootherJacobi(A, q, x, 1, 1);
solvers = [solvers; pack(jacobi, 'jacobi')];


gs = @(x) smootherGS(A, q, x, 1);
solvers = [solvers; pack(gs, 'gauss-seidel')];


mgsolver = @(x, varargin) solveMultigrid(G, A, q, 'x0', x,...
                                        'iterations', 1,...
                                        'coarsefactor', 5, ...
                                        'presmooth', 0, varargin{:});
solvers = [solvers; pack(...
           @(x) mgsolver(x, 'interpolator', 'unitary'),...
            'v-cycle unitary')];

solvers = [solvers; pack(...
           @(x) mgsolver(x, 'interpolator', 'linear'),...
            'v-cycle linear')];

[x, res, meta, partitions] = mgsolver(0*q, 'interpolator', 'mstpfa');
solvers = [solvers; pack(...
           @(x) mgsolver(x, 'interpolator', 'mstpfa', 'partitions', partitions),...
            'v-cycle MsTPFA')];



%%
n = size(A, 1) - 2;
% n = 5;
for i = 1:numel(solvers)    
    e = eigs(solvers{i}.f, size(A, 1), n);
    solvers{i}.eigen = e;
end

%%
close all
e = cellfun(@(x) x.eigen, solvers, 'UniformOutput', false);
e = horzcat(e{:});

% mean(e)
names = cellfun(@(x) x.name, solvers, 'UniformOutput', false);
plot(real(e), imag(e), '+')
legend(names)
hold on

plotCircle(1, '--k', 'linewidth', 2);
axis equal
grid on

figure;
ve = abs(e);
if 1
    % Condition number
    name = 'Operator condition number';
    dat = max(ve)./min(ve);
else
    name = 'Operator spectral radius';
    % Spectral radius
    dat = max(ve);
end
bar(dat)
set(gca, 'XTickLabel', names)
set(gca, 'yscale', 'log')
set(gcf, 'renderer', 'painters')
title(name)