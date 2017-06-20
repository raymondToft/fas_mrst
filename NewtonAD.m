function [runTime, res, nit, sol] = NewtonAD(newModel,constraints)
%% Set model
mode = 'newModel';
  initMode = struct('mode', mode, 'model',newModel);
  
% Initiate complete model
model = initiateModel(initMode);

%% Impose vertical equilibrium
[z_0, z_max] = deal(0, max(model.grid.cells.centroids(:,3)));
equil  = ode23(@(z,p) model.g .* model.oil.rhoO(p), [z_0, z_max], model.rock.p_r);
p_init = reshape(deval(equil, model.grid.cells.centroids(:,3)), [], 1);  clear equil
sW_init = zeros(model.grid.cells.num, 1);

%% Initialize for solution loop
[p_ad, sW_ad] = initVariablesADI(p_init, sW_init);

numSteps = constraints.numSteps;                 % number of time-steps
totTime  = constraints.totTime;             % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = constraints.tol;                % Newton tolerance
maxits   = constraints.maxits;                  % max number of Newton its

model.well.inRate = 1*sum(model.rock.pv(p_init))/totTime;
model.well.outRate = 0.5*model.well.inRate;

sol = repmat(struct('time',[],'pressure',[], 's', []),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', double(p_ad), ...
    's', double(sW_ad));

%% Main loop
t = 0; step = 0;
% hwb = waitbar(t,'Simulation ..');
tic;
while t < totTime
   t = t + dt;
   step = step + 1;
%    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
%       step, convertTo(t - dt, day), convertTo(t, day));
   
   [p_ad, sW_ad,nit,res] = ...
    newtonTwoPhaseAD(model,p_ad,sW_ad,p_ad,sW_ad,tol,maxits,dt);

   if nit > maxits
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, ...
                            'pressure', double(p_ad), ...
                            's', double(sW_ad));
%       waitbar(t/totTime,hwb);
   end
%      fprintf('Iteration %3d:  Res = %.4e \n', nit, res)

end
runTime = toc;
end