mrstModule add coarsegrid;

close all force;
clear;

%% Set up model

% Set up model geometry
[nx,ny,nz] = deal( 16,  16, 2);
[Dx,Dy,Dz] = deal(500, 500, 50);
grid = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
grid = computeGeometry(grid);

%plotGrid(grid); view(3); axis tight

%Set V- or F-cycle 1 == 'true';
F_cycle = 1;


% Set rock properties
  homogeneous = 'true';
  perm = 30*milli*darcy; 
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
  muO   = 5*centi*poise;
  % Compressibility range: {slighly: 10^-5 to 10^-6, compressible: 10^-3 to
  % 10^-4}psi^-1
  co      = 1e-3/barsa; %1e-4
  rho_ro = 1050*kilogram/meter^3; % 850
  rhoOS  = 750*kilogram/meter^3; % 750
    
  newModel = struct('grid',grid, 'muW',muW, 'cw',cw, 'rho_rw',rho_rw,  ...
                    'rhoWS',rhoWS, 'muO', muO, 'co',co, 'rho_ro',rho_ro, ...
                    'rhoOS', rhoOS, 'cr',cr,'p_r',p_r, 'poro',poro,'perm', perm, ...
                    'homogeneous',homogeneous, 'perm_range',perm_range, ... 
                    'gauss_filter_size',gauss_filter_size, 'std',std);
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

numSteps = 100;                 % number of time-steps
totTime  = 365*day;             % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxits   = 20;                  % max number of Newton its


model.well.inRate = 1*sum(model.rock.pv(p_init))/totTime;
model.well.outRate = 0.5*model.well.inRate;

sol = repmat(struct('time',[],'pressure',[], 's', []),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', double(p_ad), ...
    's', double(sW_ad));

%% Initiate multigrid variables

% Presmoothing steps
v1_iter = 1;

% Postsmoothing steps
v2_iter = 5;

% Number of levels

if(model.grid.cartDims(3)>1)
  k_level = floor(log(model.grid.cells.num) /log(2^3));
else
  k_level = floor(log(model.grid.cells.num) /log(2^2));
end

if(k_level < 3 || F_cycle == 1)
    % V -cycle
    cycle_index = 1;
else
    % F - cycle
    cycle_index = 2:k_level-1;
end

%% Main loop
t = 0; step = 0;
%hwb = waitbar(t,'Simulation ..');
%% Start stopwatch
tic
while t < totTime
   t = t + dt;
   step = step + 1;
%    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
%       step, convertTo(t - dt, day), convertTo(t, day));

  % Multigrid
  
%   [p_ad, sW_ad,nit] = multigridCycleV3(v1_iter,v2_iter,model,p_ad,sW_ad,tol,maxits,g,dt);
   
  [p_ad, sW_ad,nit] = FASCycle(v1_iter,v2_iter,model,p_ad,sW_ad,tol,maxits,dt,k_level,cycle_index);
   
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
     % waitbar(t/totTime,hwb);
   end
end
%% stop stopwatch

fprintf(' %f \n',toc);

%close(hwb);

%% Plot pressure evolution

% for i = 1:numSteps
%     figure(1); clf
%     subplot(2, 1, 1)
%     plotCellData(grid, sol(i).pressure);
%     colorbar (); axis equal tight ; view (3);
%     title('Pressure')
%     view(30, 40);
%     subplot(2, 1, 2)
%     plotCellData(grid, sol(i).s);
%     colorbar (); axis equal tight ; view (3);
%     caxis([0, 1])
%     view(30, 40);
%     title('Watersaturation')
%     drawnow
% end


%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
