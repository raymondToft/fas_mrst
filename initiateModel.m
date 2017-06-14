function [model] = initiateModel(initMode)
  %% Function description
  %
  % PARAMETERS:
  % G        - initialized grid structure
  % initMode - struct('mode', mode, 'model',model);
  %             mode = 'newModel': Initializes a new model with parameters from the
  %              -> model = struct(    'grid',grid 
  %                                    'muW',muW,
  %                                    'cw',cw,
  %                                    'rho_rw',rho_rw,
  %                                    'rhoWS',rhoWS 
  %                                    'muO', muO,
  %                                    'co',co,
  %                                    'rho_ro',rho_ro,
  %                                    'rhoOS', rhoOS,
  %                                    'cr',cr,
  %                                    'p_r',p_r,
  %                                    'poro',poro, 
  %                                    'perm', perm,
  %                                    'homogeneous',['true' or 'false'], 
  %                                    'permRange',permRange 
  %                                    'gaussFilterSize',gaussFilterSize,
  %                                    'std',std)
  %             mode = 'coarseModel': Reinitalization of the model
  %              -> initModel = struct( 'oldModel',oldModel, : model struct given from initial call to this
  %                                         function
  %                                      'weighting', weighting
  %
  % RETURNS:
  % p_ad_coarse - Coarsed version of the pressure eqs. with reinitialized
  %               ADI structure
  % sW_ad_coarse - Coarsed version of the pressure eqs. with reinitialized
  %               ADI structure
  % defect_p_ad_coarse - Coarsed version of the pressure defect
  % defect_s_ad_coarse - Coarsed version of the saturation defect
  % pIx_coarse   - Coarsed version of the pressure index
  % sIx_coarse   - Coarsed version of the saturation index
  %
  % COMMENTS:
  % - This coarsening function is currently not optimized for performance
  % - The function may be bugged in current state
  %
  % SEE ALSO:
  %

  
  model = 0;
  if(strcmp(initMode.mode,'newModel'))
      model = initNewModel(initMode.model);
  else
      model = initCoarseModel(initMode.model,initMode.weighting);
  end
  
end

function [model] = initNewModel(model)
  %% Define wells
  injIndex = 1;
  prodIndex = model.grid.cells.num;

  inRate = 1;
  outRate = 0.5;
  
  well = struct('injIndex',injIndex, 'prodIndex',prodIndex, 'inRate', inRate, 'outRate', outRate);
  
  %% Rock model  
  if(strcmp(model.homogeneous, 'true'))
    rock = makeRock(model.grid, model.perm, model.poro);

  elseif(strcmp(model.homogeneous, 'spe10'))
      rock = struct('poro',model.poro,'perm', model.perm);
  else
        p = gaussianField(model.grid.cartDims, model.perm_range, model.gauss_filter_size, model.std);
        K = p.^3. * (1e-5)^2./(0.81 * 72 * (1-p).^2);
        rock = makeRock(model.grid, K(:), p (:));
        figure;
        plotCellData (model.grid , convertTo ( rock.perm , milli * darcy ));
        colorbar (); axis equal tight ; view (3);
  end
  
  % Compressibility
  cr   = model.cr;
  % Reference pressure
  p_r  = model.p_r;
  pv_r = poreVolume(model.grid, rock);
  pv   = @(p) pv_r .* exp( cr * (p - p_r) );

  
  % redefine rock struct
  rock = struct('perm',rock.perm,'poro',rock.poro, ...
      'cr', cr, 'p_r',p_r, 'pv_r', pv_r, 'pv',pv);
  %% Define model for two-phase compressible fluid
  % Define a water phase
  muW    = model.muW;
  cw     = model.cw;
  rho_rw = model.rho_rw;
  rhoWS  = model.rhoWS;
  rhoW   = @(p) rho_rw .* exp( cw * (p - p_r) );
  krW = @(S) S.^2;
  
  water = struct('muW', muW, 'cw', cw, 'rho_rw', rho_rw, 'rhoWS', rhoWS, 'rhoW', rhoW, 'krW', krW);
    
  % Define oil phase
  muO   = model.muO;
  co      = model.co; 
  rho_ro = model.rho_ro; 
  rhoOS  = model.rhoOS; 
  krO = @(S) S.^3;

  rhoO   = @(p) rho_ro .* exp( co * (p - p_r) );

  oil = struct('muO', muO, 'co', co, 'rho_ro', rho_ro, 'rhoOS', rhoOS, 'rhoO', rhoO, 'krO', krO);
  
  
  %% Compute transmissibilities
  N  = double(model.grid.faces.neighbors);
  intInx = all(N ~= 0, 2);
  N  = N(intInx, :);                          % Interior neighbors
  hT = computeTrans(model.grid, rock);                 % Half-transmissibilities
  cf = model.grid.cells.faces(:,1);
  nf = model.grid.faces.num;
  T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]); % Harmonic average
  T  = T(intInx);                             % Restricted to interior
  
  %% Define discrete operators
  n = size(N,1);
  C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, model.grid.cells.num);
  grad = @(x)C*x; % Discrete gradient
  div  = @(x)-C'*x; % Discrete divergence
  avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2))); % Averaging
  upw = @(flag, x) flag.*x(N(:, 1)) + ~flag.*x(N(:, 2)); % Upwinding 

  gradz  = grad(model.grid.cells.centroids(:,3));

  operator = struct('grad', grad, 'div', div, 'avg', avg, 'upw', upw, 'gradz', gradz, 'C',C);
  
  %% Remaining variables
  %Note: Write a better description
  nc = model.grid.cells.num;
  pIx = 1:nc;
  sIx = (nc+1):(2*nc);
  p_ad = 0;
  sW_ad = 0;
  gravity reset on, g = norm(gravity);
  residual = 99e10;
  %% Multigrid cycle variables
  cycle = struct('v1', model.cycle.v1, 'v2',model.cycle.v2, ...
      'type',model.cycle.type,'level',0,'index',1,'grids',1);
  
  %% Place all model parts and help function i a "modelstruct"
  model = struct('grid',model.grid,'rock', rock, 'water', water, 'oil',oil, 'T', T, ...
      'operator', operator, 'well', well, 'pIx', pIx,'sIx',sIx, ...
      'p_ad',p_ad,'sW_ad',sW_ad,'g',g, 'residual', residual, 'cycle', cycle,'nit',0);

  
end

function [coarseModel] = initCoarseModel(coarseModel,weighting)
    %% Coarsen rock
    
    coarseModel.well.prodIndex = coarseModel.grid.cells.num;
    coarseModel.rock.perm = accumarray(coarseModel.grid.partition,coarseModel.rock.perm)./weighting;
    coarseModel.rock.poro = accumarray(coarseModel.grid.partition,coarseModel.rock.poro)./weighting;
    
        
    %Set the rock at the uncorarsened wells
    coarseModel.rock.perm(coarseModel.well.injIndex) = coarseModel.rock.perm(coarseModel.well.injIndex);
    coarseModel.rock.perm(coarseModel.well.prodIndex) = coarseModel.rock.perm(coarseModel.well.prodIndex);
      
    coarseModel.rock.poro(coarseModel.well.injIndex) = coarseModel.rock.poro(coarseModel.well.injIndex);
    coarseModel.rock.poro(coarseModel.well.prodIndex) = coarseModel.rock.poro(coarseModel.well.prodIndex);
    coarseModel.rock.pv_r = poreVolume(coarseModel.grid, coarseModel.rock);
    coarseModel.rock.pv_r(coarseModel.well.injIndex) = coarseModel.rock.pv_r(coarseModel.well.injIndex);    
    coarseModel.rock.pv_r(coarseModel.well.prodIndex) = coarseModel.rock.pv_r(coarseModel.well.prodIndex);
    coarseModel.rock.pv = @(p) coarseModel.rock.pv_r .* exp( coarseModel.rock.cr * (p - coarseModel.rock.p_r) );
      
    %% Recompute transmissibilities
    N  = double(coarseModel.grid.faces.neighbors);
    intInx = all(N ~= 0, 2);
    N  = N(intInx, :);                          % Interior neighbors
    hT = computeTrans(coarseModel.grid, coarseModel.rock);                 % Half-transmissibilities
    cf = coarseModel.grid.cells.faces(:,1);
    nf = coarseModel.grid.faces.num;
    T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]); % Harmonic average
    T  = T(intInx);                             % Restricted to interior
 
    coarseModel.T = T;
    
    %% Redefine discrete operators
    n = size(N,1);
    C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, coarseModel.grid.cells.num);
    grad = @(x)C*x; % Discrete gradient
    div  = @(x)-C'*x; % Discrete divergence
    avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2))); % Averaging
    upw = @(flag, x) flag.*x(N(:, 1)) + ~flag.*x(N(:, 2)); % Upwinding 

    gradz  = grad(coarseModel.grid.cells.centroids(:,3));

    operator = struct('grad', grad, 'div', div, 'avg', avg, 'upw', upw, 'gradz', gradz, 'C',C);
    
    coarseModel.operator = operator;
    %% Remaining variables
    nc = coarseModel.grid.cells.num;
    pIx = 1:nc;
    sIx = (nc+1):(2*nc);
  
    coarseModel.pIx = pIx;
    coarseModel.sIx = sIx;

end