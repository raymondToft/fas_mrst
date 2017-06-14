function [p_ad, sW_ad,nit,resNorm] =  ...
    newtonTwoPhaseAD(model,p_ad,sW_ad,p_ad_0,sW_ad_0,tol,maxits,dt,varargin)
   %% Function description
   %
   % PARAMETERS:
   % model    - System model structure with grid, rock, phases and operator
   %            substructs
   %
   % p_ad     - ADI struct for the pressure
   % s_ad     - ADI struct for the saturation
   % tol      - error tolerance for newton
   % maxits   - The maximum number of newton iterations perfomed
   % g        - gravity constant 
   % dt       - current time step
   % p_ad_0   - Pressure ADI from previous timestep
   % sW_ad_0  - Saturation ADI from previous timestep
   % varargin - boundary condition from multigrid cycle
   %
   % RETURNS:
   % p_ad     - Approximated values of the pressure stored in ADI structure
   % s_ad     - Approximated values of the saturation stored in ADI structure
   % res      - Residual from the final newton iteration
   % nit      - Number of iterations performed
   % 
   % COMMENTS:
   % - The body of this function is simply copyed out and slightly modified from the original file twoPhaseAD.m
   %
   % SEE ALSO:
   %

%% Body

   p0 = double(p_ad_0);
   sW0 = double(sW_ad_0);
   
   nit = 0;
   resNorm = 1e99;
   old_res = resNorm;
   scaling = 0;
   pMaxUpd =0;
   sWMaxUpd = 0;
    if(strcmp(model.grid.type,'generateCoarseGrid'))
          vol = sum(model.grid.cells.volumes);%max(model.grid.cells.volumes);
          scaling = dt/vol;
    end
   
  while (resNorm > tol) && (nit < maxits) && (old_res >= resNorm)
      old_res = resNorm;
     
      [water, oil] = computePhaseFlux(model,p_ad,sW_ad,p0,sW0,dt);
        
      % Check wether a defect have been passed or not.
      if(isempty(varargin) || isempty(varargin{1}))% || nit == 1)
          [water, oil] = computeBoundaryCondition(model,p_ad,sW_ad,water,oil);      
      else
          boundaryCondition = varargin{1};
          water =  boundaryCondition.water + water;
          oil =   boundaryCondition.oil + oil;
%           
%          
%           water_val = water(model.well.prodIndex).val;
%           oil_val = oil(model.well.prodIndex).val;
%           
%           water(model.well.prodIndex) = water(model.well.prodIndex) - water_val;% ... 
%               - water(model.well.prodIndex) + p_ad(model.well.prodIndex) - p_ad(model.well.prodIndex).val ;
%       
%           oil(model.well.prodIndex) = oil(model.well.prodIndex) - oil_val;% ...
%               - oil(model.well.prodIndex) + sW_ad(model.well.prodIndex)- sW_ad(model.well.prodIndex).val ;    
      end
      
      if(strcmp(model.grid.type,'generateCoarseGrid'))
          water = water.*scaling;
          oil = oil.*scaling;
      end
      
      % Collect all equations
      eqs = {oil, water};
      % Concatenate equations and solve for update:
      eq  = cat(eqs{:});
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      upd = -(J \ res); % Newton update
      
%       %% Update variables V1
%       pUpd = upd(model.pIx);
%       pUpdCheck = pUpd > p_ad.val*0.2; 
%       pUpd(pUpdCheck) = p_ad(pUpdCheck)*0.2;
%       pMaxUpd = pMaxUpd + sum(pUpdCheck);
%    
%       sWUpd = upd(model.sIx);
%       sWUpdCheck = sWUpd > 0.2;
%       sWUpd(sWUpdCheck) = 0.2;
%       sWMaxUpd = sWMaxUpd + sum(sWUpdCheck);  
%       
%       p_ad.val  = p_ad.val  + pUpd;% upd(model.pIx);
%       sW_ad.val = sW_ad.val + sWUpd;% upd(model.sIx);
%       sW_ad.val = min(sW_ad.val, 1);
%       sW_ad.val = max(sW_ad.val, 0);

      %% Update variables V2
      dp = upd(model.pIx);
      ds = upd(model.sIx);
      % Limit pressure to relative changes of 20%
      maxRelCh = 0.2;
      biggestChange = max(abs(dp./double(p_ad)), [], 2);
      w = min(maxRelCh./biggestChange, 1);
      % Multiply by relaxation factor
      dp = dp.*min(w);
      
      % Limit saturation to relative changes of 20%
      maxRelCh = 0.2;
      biggestChange = max(abs(ds), [], 2);
      w = min(maxRelCh./biggestChange, 1);
      % Multiply by relaxation factor
      ds = ds.*min(w);

      p_ad.val   = p_ad.val + dp;
      sW_ad.val = sW_ad.val + ds;
      sW_ad.val = min(sW_ad.val, 1);
      sW_ad.val = max(sW_ad.val, 0);
      
%     figure
%     plot(1:model.G.cells.num,sW_ad.val);
%     title('Saturation')
      
      resNorm = norm(res);
      nit     = nit + 1;
%       fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
  end
  
  model.residual = resNorm;
  %display(model.cycle.index);
%   if(~isempty(varargin))
%       fprintf('%s \n',varargin{2});
%   end
     fprintf('Grid %d, Iteration %3d:  Res = %.4e \n',model.grid.cartDims(1), nit, resNorm)
%    if(pMaxUpd > 0 || sWMaxUpd > 0)
%    fprintf('  pMaxUpd: %d, sWMaxUpd: %d \n', pMaxUpd, sWMaxUpd);
%    end
%    fprintf('Iterantions: %3d\n', nit)
end
