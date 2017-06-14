 [G, A, q, rock] = getTestCase('dims', [100 1], ...
                            'grid', 'cart',...
                            'perm', 'homo', ...
                            'stencil',  'tpfa', ...
                            'flow',   'none', ...
                            'griddims',  [1 1], ...
                            'ordering', 'default');
qnorm = norm(q, 2);
if qnorm == 0
    qnorm = 1;
end
saveplot = @(name) export_fig(['../figures/smoothers/', name, '.eps'], '-transparent');
%%
figure(2); clf;
plotCellData(G, A\q, 'edgec', 'w', 'edgea', .5)
%%

getData = @(K) sin(X*pi*K);

jacobi = @(x) smootherJacobi(A, q, x, 1, 1);
jacobi_optimal = @(x) smootherJacobi(A, q, x, 1, 2/3);
jacobi_half = @(x) smootherJacobi(A, q, x, 1, 1/2);

gs = @(x) smootherGS(A, q, x, 1);


packSolver = @(fun, name) struct('f', fun, 'name', name);

problems = {};
problems = [problems; packSolver(jacobi, 'Jacobi (\omega = 1)')];
problems = [problems; packSolver(jacobi_optimal, 'Jacobi (\omega = 2/3)')];
problems = [problems; packSolver(jacobi_half, 'Jacobi (\omega = 1/2)')];

problems = [problems; packSolver(gs, 'Gauss-Seidel')];

X = G.cells.centroids(:, 1);
packData =  @(K, name) struct('data', getData(K), 'name', ['K = ', num2str(K)], 'k', K);

data = {};
data = [data; packData(25)];
% data = [data; packData(1/2)];
% for i = 1:5:floor(G.cells.num/2)
%     data = [data; packData(i)];
% end


residuals = [];
names = {};
for d = 1:numel(data)
%     x0 = rand(G.cells.num, 1);
    x0 = data{d}.data;
    dn = data{d}.name;
    
    N = 500;
    

    for i = 1:numel(problems);
        x = x0;
        f = problems{i}.f;
        pn = problems{i}.name;
        
        res = nan(N, 1);
        for j = 1:N
            res(j) = norm(A*x - q, 2)/qnorm;
            x =  f(x);
        end
        residuals = [residuals, res]; %#ok
        names = [names; [pn, ', ', dn] ];
    end
end

%%
close all
df = get(0, 'defaultfigureposition');
figure('position', df.*[1 1 1.5 .75]);

loglog(residuals, '--', 'LineWidth', 2)
legend(names)

xlabel('# Iterations');
ylabel('||Ax-q ||_2')
grid on
axis tight
saveplot(['convergence', num2str(data{1}.k)])
%%
% Different smoothing values
K = 1;

plot(X, sin(X*2*pi*K))

%%
close all
n = 5;

df = get(0, 'defaultfigureposition');

for k = [1, 5, 10, 25]
    x0 = getData(k);
    tmp = nan(G.cells.num, n);
    x = x0;
    for i = 1:n
        for j = 1:10
            x = jacobi_optimal(x);
        end
        tmp(:, i) = x;
    end
%     norm(x, 2)/G.cells.num
    norm(x0, 2)
    figure('position', df.*[1 1 1.5 .75]);
    hold on
    plot(x0, '--', 'linewidth', 2)
    plot(tmp, '-', 'linewidth', 2)
    names = arrayfun(@(x) [num2str(x),' Jacobi'], (1:n).*10, 'UniformOutput', false);
    names = ['Initial data', names];
    legend(names, 'location', 'northeast')
    grid on
    title(['K = ', num2str(k)]);
    saveplot(['simpleLinear', num2str(k)]);    
end
%%
x = (0:0.001:1) .';

low = sin(x*pi*2);
hi = sin(x*pi*8);

data = [low, hi];

subplot(2,1,1)
plot(x, data)

offset = 3;
subs = 1:offset:numel(x)-(offset+1);

subplot(2,1,2)
plot(x(subs), data(subs, :))
% subplot(2,2,3)
% plot(X(subs), low(subs))
% subplot(2,2,4)
% plot(X(subs), hi(subs))