function precond = makeMultigridPreconditioner(G, A, q, varargin)
    % Set up solver interface
    solver = @(q, varargin) solveMultigrid(G, A, q, varargin{:});
    % Do not solve, but do generate partitions. No-op if partitions is
    % already provided.
    if mod(numel(varargin), 2) == 1
        partitions = varargin{1};
        varargin = varargin(2:end);
    else
        [x, res, meta, partitions] =  solver(q, 'doSolve', false, varargin{:}); %#ok
    end
    % Return a function handle suitable for MATLABs linear solvers
    precond = @(q, transp) solver(q, 'partitions', partitions, varargin{:});
end