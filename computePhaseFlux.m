function [water, oil] = computePhaseFlux(model,p_ad,sW_ad,p0,sW0,dt)
  % Evaluate properties
  rW = model.water.rhoW(p_ad);
  rW0 = model.water.rhoW(p0);
  rO = model.oil.rhoO(p_ad);
  rO0 = model.oil.rhoO(p0);
      
  % Define pressure drop over interface for both phases
  dp = model.operator.grad(p_ad);
  dpW = dp - model.g*model.operator.avg(rW).*model.operator.gradz; % Water
  dpO = dp - model.g*model.operator.avg(rO).*model.operator.gradz; % Oil
  % Pore volume of cells at current pressure and previous pressure 
  vol0 = model.rock.pv(p0);
  vol = model.rock.pv(p_ad);
  % Mobility: Relative permeability over constant viscosity
  mobW = model.water.krW(sW_ad)./model.water.muW;
  mobO = model.oil.krO(1-sW_ad)./model.oil.muO;
  % Define phases fluxes. Density and mobility is taken upwinded (value
  % on interface is defined as taken from the cell where the phase flow
  % is currently coming from). This gives more physical solutions than
  % averaging or downwinding.
  vW = -model.operator.upw(double(dpW) <= 0, rW.*mobW).*model.T.*dpW;
  vO = -model.operator.upw(double(dpO) <= 0, rO.*mobO).*model.T.*dpO;
  % Conservation of water and oil
  water = (1/dt).*(vol.*rW.*sW_ad - vol0.*rW0.*sW0) + model.operator.div(vW);
  oil   = (1/dt).*(vol.*rO.*(1-sW_ad) - vol0.*rO0.*(1-sW0)) + model.operator.div(vO);
end