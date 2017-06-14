% Set up the polymer equations
% In this function, the modeling equations (in this case mass conservation of oil, water and polymer) are assembled. The function returns the residual vector (one value per cell) and the Jacobians with respect to the main variables (pressure, saturation and polymer concentration).

function eqs = equationOWPolymer(state0, state, dt, bc, system, param)

   sys = system; % short-cut
   f = sys.fluid;
   poro = sys.rock.poro;

% Initialize then main variables as Automatic Differentiation (AD) variables

   p    = state.pressure;
   sW   = state.s(:,1);
   c    = state.c;
   cmax = state.cmax;
   nc   = numel(p);
   [p, sW, c] = initVariablesADI(p, sW, c);

% Initialize previous time step variables.

   p0  = state0.pressure;
   sW0 = state0.s(:,1);
   c0  = state0.c;
   cmax0 = state0.cmax;

% Compute formation volume factors

   bW     = f.bW(p);
   bO     = f.bO(p);
   bc_bW  = f.bW(bc.dirichlet.pressure); % boundary values
   bc_bO  = f.bO(bc.dirichlet.pressure); %

% Computation of the pressure gradients

   dpO = sys.grad(p, bc.dirichlet.pressure);
   dpW = dpO; % no capillary pressure

% Computation of faces formation volume factor, given as average.

   fbW = sys.averageFaceValues(bW, bc_bW);
   fbO = sys.averageFaceValues(bO, bc_bO);

% Computation of the mobilities, in the cells. The mobilities at the dirichlet boundary have been already computed.

% For an implicit solver, we compute the mobility using the current values for pressure, saturation and polymer concentration. For an explicit solver, we use the previous values. Note that this is the only modifications which are needed in the code to switch between an implicit and explicit transport solver.

   if system.implicit_transport
      [mobW, mobO, mobP] = computeMobilities(p, sW, c, cmax, f);
   else
      [mobW, mobO, mobP] = computeMobilities(p0, sW0, c0, cmax0, f);
   end

% Computation of upstream phases mobilities.

   wupc = (double(dpW)>=0);
   oupc = (double(dpO)>=0);
   fmobW = sys.upwindFaceValues(wupc, mobW, bc.mobW);
   fmobO = sys.upwindFaceValues(oupc, mobO, bc.mobO);
   fmobP = sys.upwindFaceValues(wupc, mobP, bc.mobP);

% Computation of fluxes

   vW = fbW.*fmobW.*sys.T.*dpW;
   vO = fbO.*fmobO.*sys.T.*dpO;
   vP = fbW.*fmobP.*sys.T.*dpW;

% Computation of accumulation terms for oil, water and polymer The accumulation term for polymer includes adsorption

   accumO = (sys.pv/dt).*bO.*(1 - sW);
   accumW = (sys.pv/dt).*bW.*sW;
   accumP = (sys.pv/dt).*bW.*sW.*c + (sys.pv/dt).*(f.rhoR* ( (1-poro)./poro)).* ...
            f.effads(c, cmax);

% Accumulation term for previous time step.

   accumO0 = 1/dt*bO.*sys.pv.*(1 - sW0);
   accumW0 = 1/dt*bW.*sys.pv.*sW0;
   accumP0 = (sys.pv/dt).*bW.*sW0.*c0 + (sys.pv/dt).*(f.rhoR* ( (1-poro)./poro)).* ...
            f.effads(c0, cmax0);

%Compute source term.

   sInj = bc.injection.s;
   cInj = bc.injection.c;
   influxO = zeros(nc, 1);
   influxW = zeros(nc, 1);
   influxP = zeros(nc, 1);
   influxO(bc.injection.cells) = (1 - sInj)*bc.injection.rate;
   influxW(bc.injection.cells) = sInj*bc.injection.rate;
   influxP(bc.injection.cells) = cInj*sInj*bc.injection.rate;

%assemble the equations as the sum of accumulation terms, flux terms and influxes.

   eqs{1} = accumO - accumO0 + sys.div(vO) - influxO;
   eqs{2} = accumW - accumW0 + sys.div(vW) - influxW;
   eqs{3} = accumP - accumP0 + sys.div(vP) - influxP;

end

