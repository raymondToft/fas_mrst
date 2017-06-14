mrstModule add coarsegrid;

% close all;

%% Set up model
 % Set up model geometry
[nx,ny,nz] = deal( 48,  48, 4);
[Dx,Dy,Dz] = deal(500, 500, 50);
grid = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
grid = computeGeometry(grid);

% plotGrid(grid); view(3); axis tight

% Set rock properties
  homogeneous = 'true';
  % (Outdated) value ranges
  % permeability range: {poor: 1-15, moderate: 15-20, good: 50-250, very
  % good: 250-1000
  % porosity range: {fair: 0.25, very low: 0.1}
  perm = 10*milli*darcy; 
  poro = 0.25;

  perm_range = [0.1 0.4];
  gauss_filter_size = [3 3 3];
  std = 2.5;
  % Compressibility: normally in the range of 10^-6 to 10^-7, assumed to be
  % constant
  cr   = 1e-6/barsa;
  % Reference pressure
  p_r  = 200*barsa;
 
% Set fluid properties
  % Water
  muW    = 1*centi*poise;
  cw     = 1e-5/barsa;
  rho_rw = 960*kilogram/meter^3;
  rhoWS  = 1000*kilogram/meter^3;
  
  % Define a lighter, more viscous oil phase with different relative
  % permeability function
  muO   = 1*centi*poise;
  % Compressibility range: {slighly: 10^-5 to 10^-6, compressible: 10^-3 to
  % 10^-4}psi^-1
  co      = 1e-3/barsa; %1e-4
  rho_ro = 1050*kilogram/meter^3; % 850
  rhoOS  = 750*kilogram/meter^3; % 750
     v1_iter = 0;              %Number of presmoothing steps
  v2_iter = 0;            %Number of postmoothing steps
  cycle_type = 'V_cycle';  % F_cycle or V_cycle
  cycle = struct('v1', v1_iter, 'v2',v2_iter,'type',cycle_type); 
 
  newModel = struct('grid',grid, 'muW',muW, 'cw',cw, 'rho_rw',rho_rw,  ...
                    'rhoWS',rhoWS, 'muO', muO, 'co',co, 'rho_ro',rho_ro, ...
                    'rhoOS', rhoOS, 'cr',cr,'p_r',p_r, 'poro',poro,'perm', perm, ...
                    'homogeneous',homogeneous, 'perm_range',perm_range, ... 
                    'gauss_filter_size',gauss_filter_size, 'std',std,'cycle',cycle);
  mode = 'newModel';
  initMode = struct('mode', mode, 'model',newModel);
  
% Initiate complete model
model = initiateModel(initMode);

% 
% p = linspace(100*barsa,220*barsa,50)';
% s = linspace(0,1,50)';
% plot(p/barsa, model.rock.pv_r(1).*exp(model.rock.cr*(p-model.rock.p_r)),'LineWidth',2);


%% Plot model for two-phase compressible fluid
% Water phase and a a lighter, more viscous oil phase with different relative
% permeability function

% figure;
% plot(p/barsa, [model.water.rhoW(p), model.oil.rhoO(p)],'LineWidth',2);
% legend('Water density', 'Oil density')
% 
% figure;
% plot(p/barsa, [model.water.krW(s), model.oil.krO(s)],'LineWidth',2);
% legend('krW', 'krO')
% 
% spy(model.operator.C)

%% Impose vertical equilibrium
[z_0, z_max] = deal(0, max(model.grid.cells.centroids(:,3)));
equil  = ode23(@(z,p) model.g .* model.oil.rhoO(p), [z_0, z_max], model.rock.p_r);
p_init = reshape(deval(equil, model.grid.cells.centroids(:,3)), [], 1);  clear equil
sW_init = zeros(model.grid.cells.num, 1);

%% Initialize for solution loop
[p_ad, sW_ad] = initVariablesADI(p_init, sW_init);

numSteps = 100;                 % number of time-steps
totTime  = 40*365*day;             % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxits   = 100;                  % max number of Newton its


model.well.inRate = 1*sum(model.rock.pv(p_init))/totTime;
model.well.outRate = 0.5*model.well.inRate;

sol = repmat(struct('time',[],'pressure',[], 's', []),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', double(p_ad), ...
    's', double(sW_ad));

%% Initiate multigrid variables


%% Main loop
t = 0; step = 0;
% hwb = waitbar(t,'Simulation ..');
tic;
while t < totTime
   t = t + dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));

   
   [p_ad, sW_ad,nit,resNorm] = ...
    newtonTwoPhaseAD(model,p_ad,sW_ad,p_ad,sW_ad,tol,maxits,dt);
%   if mod(step,10) == 0
%     figure
%     subplot(2, 1, 1); plot(model.grid.cells.indexMap,p_ad.val);
%     title('Pressure')
%     
%     subplot(2, 1, 2); plot(model.grid.cells.indexMap,sW_ad.val);
%     title('Saturation')
%     drawnow
%    end
   
   if nit > maxits
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, ...
                            'pressure', double(p_ad), ...
                            's', double(sW_ad));
%       waitbar(t/totTime,hwb);
   end
end

fprintf('Runtime: %.2f, Residual: %.4e, Iterations: %d \n', toc, resNorm,nit);

% close(hwb);
%% Plot pressure evolution
% 
% for i = 1:numSteps
%     figure(1); clf
%     subplot(2, 1, 1)
%     plotCellData(grid, sol(i).pressure);
%     title('Pressure')
%     view(30, 40);
%     subplot(2, 1, 2)
%     plotCellData(grid, sol(i).s);
%     caxis([0, 1])
%     view(30, 40);
%     title('Watersaturation')
%     drawnow
% end
