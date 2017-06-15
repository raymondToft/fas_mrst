function [p_approx, sW_approx,nit,resNorm] ...
    = FASCycle(model,p_ad, sW_ad,p_ad_0,sW_ad_0,tol,maxits,dt, varargin)
  %% Function description
  %
  % PARAMETERS:
  % model    - System model structure with grid, rock, phases and operator
  %            substructs
  %
  % p_ad     - ADI struct for the pressure
  % s_ad     - ADI struct for the saturation
  % tol      - error tolerance for newton at the most coarse grid
  % maxits   - The maximum number of newton iterations perfomed on the most
  %            coarsed grid
  % g        - gravity constant 
  % dt       - current time step
  % varargin - defect/boundary_condition from previous coarsening level
  %
  % RETURNS:
  % p_ad     - Approximated values of the pressure stored in ADI structure
  % s_ad     - Approximated values of the saturation stored in ADI structure
  % nit      - Number of newton iterations performed at the most coarse grid
  % COMMENTS:
  %   This is a Full Approximation Scheme, FAS, in progress
  %
  % SEE ALSO:
  %

  %% Presmoothing
  if(~isempty(varargin))
      % If varargin contains a cell, extract the struct from the first cell
      % position
     varargin = varargin{1}; 
  end
%   fprintf('Leve1 %d, Presmooth:  \n',model.cycle.level);
  if(model.residual < tol)
    v1 = model.cycle.v1;
  else
    v1 = model.cycle.v1 + (model.cycle.grids-model.cycle.level)*2;
  end
    [p_ad,sW_ad] = newtonTwoPhaseAD(model, p_ad,sW_ad,p_ad_0,sW_ad_0,tol,v1,dt, varargin);
  
      
  %% Find the defect
  [water, oil] = computePhaseFlux(model,p_ad,sW_ad,p_ad_0,sW_ad_0,dt);
  water = (-1)*water;
  oil = (-1)*oil;
  [water, oil] = computeBoundaryCondition(model,p_ad,sW_ad,water,oil,true);
  water.val = -water.val;
  oil.val = -oil.val;
  defect = struct('water',water,'oil', oil);
  
  %% Set up of coarse grid
  [coarse_model,coarse_p_ad, coarse_sW_ad, coarse_p_ad_0, coarse_sW_ad_0,coarse_defect] ...
    = coarsening(model, p_ad, sW_ad, p_ad_0, sW_ad_0, defect);
  
  %% Compute new right hand side
  [coarse_water, coarse_oil] = computePhaseFlux(coarse_model,coarse_p_ad,coarse_sW_ad,coarse_p_ad_0,coarse_sW_ad_0,dt);
  
  rhs_water = coarse_defect.water.val + coarse_water;
  rhs_oil = coarse_defect.oil.val + coarse_oil;

  boundary_condition = struct('water',rhs_water,'oil',rhs_oil);
  %% Multigrid cycle
  if(coarse_model.cycle.level > 1)
%         fprintf('Leve1 %d, send FASCycle: \n',coarse_model.cycle.level);
        coarse_model.cycle.level = coarse_model.cycle.level - 1;
      
        [coarse_approx_p, coarse_approx_sW,nit,resNorm] = FASCycle(coarse_model, ...
          coarse_p_ad,coarse_sW_ad,coarse_p_ad_0, coarse_sW_ad_0,tol,maxits,dt,boundary_condition);
      
      coarse_model.cycle.level = coarse_model.cycle.level + 1; 
%       fprintf('Leve1 %d, return FASCycle: \n',coarse_model.cycle.level);
      
      if(strcmp(model.cycle.type,'F_cycle') && ~isempty(coarse_model.cycle.index) && coarse_model.cycle.level == coarse_model.cycle.index(1))
          coarse_model.cycle.index(1) = []; %Pop current index
%           fprintf('Leve1 %d, reSend FASCycle: \n',coarse_model.cycle.level);
          coarse_model.cycle.level = coarse_model.cycle.level - 1; 
          
          [coarse_approx_p, coarse_approx_sW,nit,resNorm] = FASCycle(coarse_model, ...
              coarse_approx_p,coarse_approx_sW,coarse_p_ad_0, coarse_sW_ad_0,tol,maxits,dt, boundary_condition);
          coarse_model.cycle.level = coarse_model.cycle.level + 1; 
%           fprintf('Leve1 %d, reReturn FASCycle: \n',coarse_model.cycle.level);
          
      end
  else
    % Multigrid core: compute a approximation on the corse grid
%     fprintf('Level %d , Solve:     \n',coarse_model.cycle.level);
    [coarse_approx_p,coarse_approx_sW,nit,resNorm] ...
          = newtonTwoPhaseAD(coarse_model,coarse_p_ad,coarse_sW_ad,coarse_p_ad_0, coarse_sW_ad_0,tol,maxits,dt,boundary_condition);
    model.residual = resNorm;
  end
  %% Compute correction
  corse_correction_p =  coarse_approx_p - coarse_p_ad.val;
  corse_correction_sW = coarse_approx_sW - coarse_sW_ad.val;
  
  %% Interpolating soluton from coarsed grid and compute ccorrected approximation
  [model,fine_correction_p, fine_correction_sW] = interpolate(model,coarse_model,  ...
        corse_correction_p , corse_correction_sW);
  
  p_ad.val =   p_ad.val + fine_correction_p;
  sW_ad.val = sW_ad.val + fine_correction_sW;
  
  %% Postsmoothing
  
%   fprintf('Level %d, Postsmooth: \n',model.cycle.level);
   [p_approx,sW_approx,nit,resNorm] = newtonTwoPhaseAD(model,p_ad,sW_ad,p_ad_0,sW_ad_0,tol,model.cycle.v2*(model.cycle.grids-model.cycle.level+3),dt);
   
end
