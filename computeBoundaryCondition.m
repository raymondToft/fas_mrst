function [water, oil] = computeBoundaryCondition(model,p_ad,sW_ad,water,oil,varargin)
  
  sign = 1;
%   if(~isempty(varargin))
%       % Computation of defect
%      sign = -1;
% %   elseif(strcmp(model.G.type{1},'generateCoarseGrid'))   
% %       scale = norm(model.G.cartDims)/norm(model.G.parent.cartDims); 
% %       sign = scale;
%   end

% %% Simple border condition
%      water(model.well.injIndex) = water(model.well.injIndex) - model.well.inRate.*model.water.rhoWS;
%       % Set production cells to fixed pressure of 200 bar and zero water
%       water(model.well.prodIndex) = sW_ad(model.well.prodIndex);
%       oil(model.well.prodIndex) = p_ad(model.well.prodIndex) - 200*barsa;

 %% Set injection well boundary conditions
  % Insert volumetric source term multiplied by density
%   water(model.well.injIndex) = water(model.well.injIndex) - sign*model.well.inRate.*model.water.rhoWS;
%   
%   %% Set production well boundary conditions
%   
%   rW = model.water.rhoW(p_ad);
% %   rO = model.oil.rhoO(p_ad);
% % dp = model.operator.grad(p_ad);
%   mobW = model.water.krW(sW_ad)./model.water.muW;
% %   mobO = model.oil.krO(1-sW_ad)./model.oil.muO;
%       
%   q_w = (rW(model.well.prodIndex)*model.rock.perm(model.well.prodIndex)*mobW(model.well.prodIndex)*model.grid.cells.centroids(model.well.prodIndex,3)) ...
%     *(200*barsa - p_ad(model.well.prodIndex).val);
% %   q_o = (rO(model.well.prodIndex)*model.rock.perm(model.well.prodIndex)*mobO(model.well.prodIndex)*model.grid.cells.centroids(model.well.prodIndex,3)) ...
% %     *(200*barsa - p_ad(model.well.prodIndex).val);
%   
%   % The ADI jacobi variables are also set to the production cells
%   water(model.well.prodIndex) = water(model.well.prodIndex) - water(model.well.prodIndex) ...
%       + p_ad(model.well.prodIndex) - p_ad(model.well.prodIndex).val - sign*q_w.val;
%   oil(model.well.prodIndex) = oil(model.well.prodIndex) - oil(model.well.prodIndex) ...
%       + sW_ad(model.well.prodIndex)- sW_ad(model.well.prodIndex).val;% - sign*q_o.val;
%    
  
  %% Default boundary condition
  % Insert volumetric source term multiplied by density
  water(model.well.injIndex) = water(model.well.injIndex) - model.well.inRate.*model.water.rhoWS;
  % Set production cells to fixed pressure of 200 bar and zero water
  water(model.well.prodIndex) = sW_ad(model.well.prodIndex);
  oil(model.well.prodIndex) = p_ad(model.well.prodIndex) - 200*barsa;

end