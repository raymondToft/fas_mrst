% compute the mobilities
% This function computes the oil, water and polymer mobilities given pressure, saturation and concentration values.

function [mobW, mobO, mobP] = computeMobilities(p, sW, c, cmax, fluid)
% Setup some polymer parameters
% Mixing coefficient ? is stored in mixpar

   mixpar = fluid.mixPar;

   cbar   = c/fluid.cmax;

% The variable muM denote the viscosity multplication factor, that is, ?m(c)/?w see tutorial polymer note.

   muM = @(conc) (fluid.muWMult(conc));

% Maximum values for the viscosity multipliers.

   muMmax = muM(fluid.cmax);
   muMmin = muM(0);

% Computing effective viscosity multiplicators for both phases, which correspond to ?w,e/?w and ?p,eff/?w in tutorial polymer note.

   conc_mult = muM(c).^(mixpar).*(muMmax^(mixpar - 1)*cbar + muMmin^(mixpar - 1)*(1 - ...
                                                     cbar)).^ (-1);
   concP_mult = (cbar + (muMmax/muMmin)^(1 - mixpar)*(1-cbar)).^(-1);

% Add permeability reduction effect, see Rk.

   permRed = 1 + ((fluid.rrf-1)/fluid.adsMax)*fluid.effads(c, cmax);
   conc_mult  = conc_mult.*permRed;

% Compute the mobilities for oil, water and polymer

   [krW, krO] = fluid.relPerm(sW);
   muW = fluid.muW(p);
   muW = conc_mult.*muW;
   mobW = krW./muW;

   mobP = c.*concP_mult.*mobW;

   muO = fluid.BOxmuO(p)./fluid.BO(p);
   mobO = krO./muO;

end

