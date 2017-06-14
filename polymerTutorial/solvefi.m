%%
%
% perlOneLineDescription(Single step non-linear solver)
% 
% Compute the next |state| of the system for given previous state |state0| and time step
% |dt|.
%
% The discretized residual equation are assembled in the function
% perlAddLink(equationOWPolymer).
%

function [state, convergence, its] = solvefi(state0, dt, bc, system, varargin)

   opt = struct('verbose', false);
   opt = merge_options(opt, varargin{:});

   fluid = system.fluid;
   
   meta = struct('converged'  , false                       , ...
                 'stopped'    , false                       , ...
                 'relax'      , 1                           , ...
                 'stagnate'   , false                       , ...
                 'iteration'  , 0                           , ... 
                 'res_history', []                          );

   timer = tic;
   
   converged = false;
   state = state0;

   if mrstVerbose
      fprintf('%13s%-26s%-36s\n', '', 'CNV (oil, water)', 'MB (oil, water)');
   end

   %%
   % The function perlAddLink(equationOWPolymer) assembles the equation and returns the
   % residuals and the Jacobians.
   %
   
   equation = @(state) equationOWPolymer(state0, state, dt, bc, system);
   
   %%
   % We start with the Newton iterations
   
   while ~ (converged || meta.stopped),
      
      % Save iteration number in meta info
      meta.iteration = meta.iteration + 1;

      %%
      % We assemble the equations.
      % The variable |eqs| is an Automatic Differentiation (AD) variable
      %
      
      eqs = equation(state);

      %%
      % We call a standard linear solve to compute the Newton step |dx|.
      %
      
      dx = SolveEqsADI(eqs, []);
      
      %%
      % We update |state|.
      %
      
      state      = updateState(state, dx, system);
      
      %%
      % We compute the residual values by calling |getResiduals|.
      % This function detects oscillation and stagnation
      %
      
      [meta, residuals] = getResiduals(meta, eqs, system, false);

      %%
      % We test for convergence. Here we use a very stringent test given by a max norm
      %
      
      if mrstVerbose
         fprintf('Newton iteration %d, max residual: %5g\n', meta.iteration, ...
                 max(abs(residuals)));
      end
      converged = all(residuals <= system.nonlinear.tol); 

      if meta.stagnate,
         warning('newt:stagnate', 'Non-linear solver stagnated...')
      end
      
      if meta.iteration > system.nonlinear.maxIterations
         meta.stopped = true;
      end
      
   end

   if meta.stopped,
      warning('newt:maxit', ...
              ['Non-linear solver did not converge, stopped ', ...
               'by max iterations...']);
      convergence = false;
   else
      its = meta.iteration;
      convergence = true;
   end
   
   dispif(mrstVerbose, 'Completed %d iterations in %1.2f s\n', ...
          meta.iteration, toc(timer));

end



%% 
% We define the |updateState| function which, given an increment, returns the updated state.
%

function state = updateState(state, dx, system)

   fluid = system.fluid;
   
   dp = dx{1};
   ds = dx{2};
   dc = dx{3};
   
   state.pressure = state.pressure + dp;

   %%
   % We cap the saturations and concentrations values.
   %
   
   state.s(:, 1) = max(0, min((state.s(:, 1) + ds), 1));
   state.s(:, 2) = 1 - state.s(:, 1);
   state.c = max(0, min((state.c + dc), fluid.cmax));
   state.cmax = max(state.c, state.cmax);
   
end
