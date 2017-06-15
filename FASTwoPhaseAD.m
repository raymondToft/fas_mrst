function [runTime, res, nit, sol] = FASTwoPhaseAD(newModel,constraints)
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

%% Initiate multigrid variables

% Number of levels
if(model.grid.cartDims(3)>1)
  model.cycle.level = floor(log(model.grid.cells.num) /log(2^3));
else
    model.cycle.level = floor(log(model.grid.cells.num) /log(2^2));
end
model.cycle.grids = model.cycle.level;

if(~strcmp(model.cycle.type, 'V_cycle') && ~strcmp(model.cycle.type, 'F_cycle'))
    fprintf('Error: %s cycle not supported. Please write V_cycle or F_cycle\n',model.cycle.type);
end

%Set cycle index for F-cycle. V-cycle is set as default
if(model.cycle.level < 3 || strcmp(model.cycle.type, 'F_cycle'))
    model.cycle.index = 2:model.cycle.level-1;
end

v1 = model.cycle.v1;
%% Main loop
t = 0; step = 0;
% hwb = waitbar(t,'Simulation ..');
tic;
while t < totTime
   t = t + dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));
  p_ad_0 = p_ad;
  sW_ad_0 = sW_ad;
  res = 99e10;
%   if (step==1)
%         [p_ad, sW_ad,nit,res] = newtonTwoPhaseAD(model,p_ad,sW_ad,p_ad_0,sW_ad_0,tol,maxits,dt);
%   end
%  [p_ad, sW_ad,nit,res] = newtonTwoPhaseAD(model,p_ad,sW_ad,p_ad_0,sW_ad_0,tol,10,dt);
  
%   while(model.residual > tol)
  [p_ad, sW_ad,nit,res] = FASCycle(model,p_ad,sW_ad,p_ad_0,sW_ad_0,tol,maxits,dt);
%   model.residual = res;
%   end
  if(res > tol)
    [p_ad, sW_ad,nit,res] = ...
    newtonTwoPhaseAD(model,p_ad,sW_ad,p_ad_0,sW_ad_0,tol,10+0*model.cycle.v2,dt);
  end
%   if(nit == model.cycle.v2)
%       model.cycle.v1 = v1+2;
%   else
%     model.cycle.v1 = v1;
%   end
   if nit > maxits
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, ...
                            'pressure', double(p_ad), ...
                            's', double(sW_ad));
%       waitbar(t/totTime,hwb);
   end
%    fprintf('Iteration %3d:  Res = %.4e \n', nit, res)

end
runTime = toc;
% close(hwb);

end