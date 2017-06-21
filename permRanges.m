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
  homogeneous = 'spe10';
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

diary '8_TestCase_permRanges_banach03.txt'
fprintf('This is the testing results from running FAS and Newton with different heterogeneous permeability fields');
% 
 %Test A 
 newModel.perm_range = [0.15 0.4];
 
         p = gaussianField(grid.cartDims, perm_range, gauss_filter_size, std);
         K = p.^3. * (1e-5)^2./(0.81 * 72 * (1-p).^2);
         rock = makeRock(grid, K(:), p (:));
         figure;
         plotCellData (grid , convertTo ( rock.perm , milli * darcy ));
         colorbar (); axis equal tight ; view (3);
 newModel.perm = rock.perm;
 newModel.poro = rock.poro;
         
 fprintf('\n FAS Test A: [10 - 300] \n');
 fprintf(' RunTime, Residual, Iterations \n');
 result_A_FAS = runSumulationFAS(newModel, constraints, nSimulations);
% 
% 
 fprintf('\n Newton Test A: [10 - 300] \n');
 fprintf(' RunTime, Residual, Iterations \n');
 result_A_N = runSumulationNewton(newModel, constraints, nSimulations);
 
 %Test B 
 newModel.perm_range = [0.3 0.45];
         p = gaussianField(grid.cartDims, perm_range, gauss_filter_size, std);
         K = p.^3. * (1e-5)^2./(0.81 * 72 * (1-p).^2);
         rock = makeRock(grid, K(:), p (:));
         figure;
         plotCellData (grid , convertTo ( rock.perm , milli * darcy ));
         colorbar (); axis equal tight ; view (3);
 newModel.perm = rock.perm;
 newModel.poro = rock.poro;
 
 
 fprintf('\n FAS Test B: [100 500]\n');
 fprintf(' RunTime, Residual, Iterations \n');
 result_B = runSumulationFAS(newModel, constraints, nSimulations);
 
 fprintf('\n Newton Test B: [100 500]\n');
 fprintf(' RunTime, Residual, Iterations \n');
 result_B_N = runSumulationNewton(newModel, constraints, nSimulations);
 % 
 %Test C 
 newModel.perm_range = [0.15 0.3];
         p = gaussianField(grid.cartDims, perm_range, gauss_filter_size, std);
         K = p.^3. * (1e-5)^2./(0.81 * 72 * (1-p).^2);
         rock = makeRock(grid, K(:), p (:));
         figure;
         plotCellData (grid , convertTo ( rock.perm , milli * darcy ));
         colorbar (); axis equal tight ; view (3);
 newModel.perm = rock.perm;
 newModel.poro = rock.poro;
         
 fprintf('\n FAS Test C: [10 100]\n');
 fprintf(' RunTime, Residual, Iterations \n');
 result_C = runSumulationFAS(newModel, constraints,nSimulations);
 
 fprintf('\n NewtonTest C: [10 100]\n');
 fprintf(' RunTime, Residual, Iterations \n');
 result_C_N = runSumulationNewton(newModel, constraints,nSimulations);
%Test D 
ewModel.perm_range = [0.25 0.52];
       p = gaussianField(grid.cartDims, perm_range, gauss_filter_size, std);
       K = p.^3. * (1e-5)^2./(0.81 * 72 * (1-p).^2);
       rock = makeRock(grid, K(:), p (:));
       figure;
       plotCellData (grid , convertTo ( rock.perm , milli * darcy ));
       colorbar (); axis equal tight ; view (3);
newModel.perm = rock.perm;
newModel.poro = rock.poro;
        
fprintf('\n FAS Test D: [50 1000]\n');
fprintf(' RunTime, Residual, Iterations \n');
result_D = runSumulationFAS(newModel, constraints,nSimulations);

fprintf('\n Newton Test D: [50 1000]\n');
fprintf(' RunTime, Residual, Iterations \n');
result_D_N = runSumulationNewton(newModel, constraints,nSimulations);


%Test E 
newModel.perm_range = [0.4 0.52];
        p = gaussianField(grid.cartDims, perm_range, gauss_filter_size, std);
        K = p.^3. * (1e-5)^2./(0.81 * 72 * (1-p).^2);
        rock = makeRock(grid, K(:), p (:));
        figure;
        plotCellData (grid , convertTo ( rock.perm , milli * darcy ));
        colorbar (); axis equal tight ; view (3);
 newModel.perm = rock.perm;
newModel.poro = rock.poro;

fprintf('\n FAS Test E: [300 1000]\n');
fprintf(' RunTime, Residual, Iterations \n');
result_E = runSumulationFAS(newModel, constraints, nSimulations);

fprintf('\n Newton Test E: [300 1000]\n');
fprintf(' RunTime, Residual, Iterations \n');
result_E_ = runSumulationNewton(newModel, constraints, nSimulations);

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
function [result] = runSumulationNewton(newModel,constraints,nSimulations)
    result = struct('variable',[],'runTime', [], 'res', [], 'nit',[]);
    for i = 1:nSimulations
%         result(i,:) = FASTwoPhaseAD(newModel,constraints);
        [runTime, res, nit]  = NewtonAD(newModel,constraints);
        result(i) = struct('variable',newModel.muO,'runTime', runTime, 'res', res, 'nit',nit);

        fprintf('%.2f %.4e %d \n', runTime, res, nit);
    end
end


function [result] = runSumulationFAS(newModel,constraints,nSimulations)
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
