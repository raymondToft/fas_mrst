function [coarse_model,p_ad_coarse, sW_ad_coarse, p_ad_0_coarse, sW_ad_0_coarse, coarse_defect] ...
    = coarsening(model, p_ad, sW_ad, p_ad_0, sW_ad_0, varargin)
  %% Function description
  %
  % PARAMETERS:
  % model    - System model structure with grid, rock, phases and operator
  %            substructs
  % p_ad     - ADI struct for the pressure
  % s_ad     - ADI struct for the saturation
  % defect   - The defect of the current approximization. (NOTE: varargin)
  %
  % RETURNS:
  % p_ad_coarse - Coarsed version of the pressure eqs. with reinitialized
  %               ADI structure
  % sW_ad_coarse - Coarsed version of the pressure eqs. with reinitialized
  %               ADI structure
  % defect_p_ad_coarse - Coarsed version of the pressure defect
  % defect_s_ad_coarse - Coarsed version of the saturation defect
  %
  % COMMENTS:
  % - This coarsening function is currently not optimized for performance
  % - The function may be bugged in current state
  %
  % SEE ALSO:
  %

  %% Create new coarse model
  % Set up coarse grid
  coarse_dims = ceil(model.grid.cartDims/2);
  coarse_dims(3) = model.grid.cartDims(3);
  partition  = partitionCartGrid(model.grid.cartDims,coarse_dims);
  
  if(size(partition,1) < model.grid.cells.num)
      diff = model.grid.cells.num - size(partition,1);
      partition = cat(1,partition,partition(end-diff+1:end,end));
  end
  if(partition(2) == 1 )
    partition(2:end) = partition(2:end) +1;
  end
  partition(model.well.prodIndex) = partition(model.well.prodIndex) + 1;
  
  coarse_grid = generateCoarseGrid(model.grid, partition);
  coarse_grid = coarsenGeometry(coarse_grid);
  
  % Add fields to the coarse grid to ensure that it passes as a
  % regular grid for our purposes.
  coarse_grid.cartDims = coarse_dims;
  
  weighting = accumarray(partition,1);
  
   model.grid = coarse_grid;
   mode = 'coarseModel';
   coarseModel = struct('mode', mode,'model',model,'weighting',weighting);
   coarse_model = initiateModel(coarseModel);
   
  %% Restrict AD variables and defect
  
  coarse_p_init = accumarray(partition, p_ad.val)./weighting;
  coarse_sW_init = accumarray(partition,sW_ad.val)./weighting;
  
  % Until a better aproach is found, the ADI varaables is re-initiated 
  [p_ad_coarse, sW_ad_coarse] = initVariablesADI(coarse_p_init, coarse_sW_init);

  coarse_p_0 = accumarray(partition, p_ad_0.val)./weighting;
  coarse_sW_0 = accumarray(partition,sW_ad_0.val)./weighting;
  
  % Until a better aproach is found, the ADI varaables is re-initiated 
  [p_ad_0_coarse, sW_ad_0_coarse] = initVariablesADI(coarse_p_0, coarse_sW_0);
  
  coarse_defect = 0;
  if(~isempty(varargin))
    defect = varargin{1};
    [coarse_water_defect,coarse_oil_defect] ...
      = initVariablesADI(accumarray(partition,defect.water.val),accumarray(partition,defect.oil.val));

    coarse_defect = struct('water',coarse_water_defect, 'oil',coarse_oil_defect);
  end
  
  % Well conditions
  
  coarse_model.well.inRate = model.well.inRate;
  coarse_model.well.outRate = model.well.outRate;

%   plotCoarseGrid = struct('cells',coarse_model.grid.cells,'faces',coarse_model.grid.faces,'cartDims',coarse_model.grid.cartDims,'griddim',coarse_model.grid.griddim,'type',{'tensorGrid'  'cartGrid'  'computeGeometry'});
%   plotGrid(plotCoarseGrid); view(3);

end