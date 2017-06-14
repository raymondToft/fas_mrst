mrstModule add mrst-gui coarsegrid msfvm deckformat coarsegrid mex
assert(exist('jacobiUpdateBasis', 'file') == 2)
solvers = {[], 'gmres'};
% solvers = {[]};
% solvers = {'gmres'};
for k = 1:numel(solvers)
%     for i = 1:2
%     for i = 3:4
    for i = 6
%         amgGridTesting('direct', i, 0, solvers{k});
%         disp('direct - ok')
        amgGridTesting('standard', i, 0, solvers{k});
        disp('standard - ok')
%         amgGridTesting('constant', i, 0, solvers{k});
%         disp('constant - ok')
%         amgGridTesting('constant', i, 1, solvers{k});
%         disp('constant-1 - ok')

    end
end
%%
for k = 1:numel(solvers)
    for i = 6
%         amgGridTesting('direct', i, 0, solvers{k});
%         disp('direct - ok')
%         amgGridTesting('standard', i, 0, solvers{k});
%         disp('standard - ok')
        amgGridTesting('constant', i, 0, solvers{k});
        disp('constant - ok')
        amgGridTesting('constant', i, 1, solvers{k});
        disp('constant-1 - ok')
    end
end

%%
for i = 1:6
    amgGridTesting('constant', i, 0, []);
end
