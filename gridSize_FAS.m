%% This script run tests for sensitivity of chnages in the mobility ratio
% Set up model 
% This is the default values for all simulations
% Spesific values for the test case are specifyed below with the function
% call to the FAS script

mrstModule add coarsegrid;

close all;

%% Set up model
 % Set up model geometry
[nx,ny,nz] = deal( 48,  48, 4);
[Dx,Dy,Dz] = deal(500, 500, 50);
grid = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
grid = computeGeometry(grid);

%plotGrid(grid); view(3); axis tight

% Set rock properties
  homogeneous = 'false';
  % (Outdated) value ranges
  % permeability range: {poor: 1-15, moderate: 15-20, good: 50-250, very
  % good: 250-1000
  % porosity range: {fair: 0.25, very low: 0.1}
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
  
  v1_iter = 1;              %Number of presmoothing steps
  v2_iter = 5;            %Number of postmoothing steps
  cycle_type = 'V_cycle';  % F_cycle or V_cycle
  cycle = struct('v1', v1_iter, 'v2',v2_iter,'type',cycle_type); 
  
newModel = struct('grid',grid, 'muW',muW, 'cw',cw, 'rho_rw',rho_rw,  ...
                    'rhoWS',rhoWS, 'muO', muO, 'co',co, 'rho_ro',rho_ro, ...
                    'rhoOS', rhoOS, 'cr',cr,'p_r',p_r, 'poro',poro,'perm', perm, ...
                    'homogeneous',homogeneous, 'perm_range',perm_range, ... 
                    'gauss_filter_size',gauss_filter_size, 'std',std,'cycle',cycle);

numSteps = 100;                 % number of time-steps
totTime  = 10*365*day;             % total simulation time
tol      = 1e-5;                % Newton tolerance
maxits   = 100;                  % max number of Newton its

constraints = struct('numSteps',numSteps,'totTime',totTime, 'tol',tol, 'maxits', maxits);
  
%% Run Simulations
nSimulations = 5;

diary 'TestCase_homogeneous_FAS.txt'
fprintf('This is the testing results from running FAS with different grid sizes for a heterogeneous permeability field (perm=30)');

%Test A 
fprintf('\n FAS Test A: 16*16*4 \n');
fprintf(' RunTime, Residual, Iterations \n');
[nx,ny,nz] = deal( 16,  16, 4);
grid = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
grid = computeGeometry(grid);
newModel.grid = grid;
result_A = runSumulation(newModel, constraints, nSimulations);

%Test B 
fprintf('\n FAS Test B: 24*24*4\n');
fprintf(' RunTime, Residual, Iterations \n');
[nx,ny,nz] = deal( 24,  24, 4);
grid = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
grid = computeGeometry(grid);
newModel.grid = grid;
result_B = runSumulation(newModel, constraints, nSimulations);

%Test C 
fprintf('\n FAS Test C: 32*32*4\n');
fprintf(' RunTime, Residual, Iterations \n');
[nx,ny,nz] = deal( 32,  32, 4);
grid = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
grid = computeGeometry(grid);
newModel.grid = grid;
result_C = runSumulation(newModel, constraints,nSimulations);

%Test D 
fprintf('\n FAS Test D: 48*48*4\n');
fprintf(' RunTime, Residual, Iterations \n');
[nx,ny,nz] = deal( 48,  48, 4);
grid = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
grid = computeGeometry(grid);
newModel.grid = grid;
result_D = runSumulation(newModel, constraints,nSimulations);

%Test E 
fprintf('\n FAS Test E: 64*64*4\n');
fprintf(' RunTime, Residual, Iterations \n');
[nx,ny,nz] = deal( 64,  64, 4);
grid = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
grid = computeGeometry(grid);
newModel.grid = grid;
result_E = runSumulation(newModel, constraints, nSimulations);


% 
% fprintf('MobilityRatio Maximum Average relAverage Variance relVariance')
% 
% refineAndPrintResults(result_A,'GridSize16')
% 
% refineAndPrintResults(result_B,'MobilityRatio')
% 
% refineAndPrintResults(result_C,'MobilityRatio')
% 
% refineAndPrintResults(result_D,'MobilityRatio')
% 
% refineAndPrintResults(result_E,'MobilityRatio')
% 
% refineAndPrintResults(result_F,'MobilityRatio')

diary off


figure; % Lazy way of getting notice of the end of the simulations
function [result] = runSumulation(newModel,constraints,nSimulations)
    result = struct('variable',[],'runTime', [], 'res', [], 'nit',[]);
    for i = 1:nSimulations
%         result(i,:) = FASTwoPhaseAD(newModel,constraints);
        [runTime, res, nit]  = FASTwoPhaseAD(newModel,constraints);
        result(i) = struct('variable',newModel.muO,'runTime', runTime, 'res', res, 'nit',nit);

        fprintf('%.2f %.4e %d \n', runTime, res, nit);
    end
end

function [] = refineAndPrintResults(results, variable)
    average = mean(results.runTime);
    maximum = max(results.runTime);
    variance = var(results.runTime); 
    relAverage = average/maximum;
    relVariance = variance/maximum;
    fprintf('%s: %d %f %f %f %f %f \n',variable, results.variable,maximum, average, relAverage, ...
    variance, relVariance);
end
