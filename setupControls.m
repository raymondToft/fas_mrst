% Set up the controls
% There are two types of controls:
%       Dirichlet boundary conditions, where the pressure, saturation and polymer concentration are given.
%       Injection fluxes, where sources are included at given cells. We have thus not included well models in this tutorial.
% The controls locations have been hard coded here, but they can be changed.
% We have implemented controls for two simulation cases 1D and 2D

function bc = setupControls(G, bc, sim_case)

   switch sim_case

% We start with the 1D case

     case '1D'

% We consider given pressure, saturation and concentration on the outer right face

       dirichlet_cell = G.cells.num; % Corresponds cell number of outer right cell for the
                                     % Cartesian grid we have defined.
       dirichlet_faces = G.cells.faces(G.cells.facePos(dirichlet_cell): ...
                                       G.cells.facePos(dirichlet_cell + 1) - 1 , :);

       bc.dirichlet.faces =  dirichlet_faces((dirichlet_faces(:, 2) == 2), 1);
       ndf = numel(bc.dirichlet.faces);
       bc.dirichlet.pressure = bc.dirichlet.pressure*ones(ndf, 1);
       bc.dirichlet.s = bc.dirichlet.s*ones(ndf, 1);
       bc.dirichlet.c = bc.dirichlet.c*ones(ndf, 1);

% We consider given influx on the outer left cell.

       injection_cells = 1; % Corresponds to cell number of outer left cell for the
                            % Cartesian grid we have defined.
       bc.injection.cells = injection_cells;

% We set up the 2D case

     case '2D'

% We consider given pressure, saturation and concentration on the faces of the upper-right cell.

       dirichlet_cell = G.cells.num; % Corresponds to upper-right cell for the Cartesian grid we
                                     % have defined.
       dirichlet_faces = G.cells.faces(G.cells.facePos(dirichlet_cell): ...
                                       G.cells.facePos(dirichlet_cell + 1) - 1 , :);

       bc.dirichlet.faces =  dirichlet_faces((dirichlet_faces(:, 2) == 2) | (dirichlet_faces(:, 2) == ...
                                                         4) , 1);
       ndf = numel(bc.dirichlet.faces);
       bc.dirichlet.pressure = bc.dirichlet.pressure*ones(ndf, 1);
       bc.dirichlet.s = bc.dirichlet.s*ones(ndf, 1);
       bc.dirichlet.c = bc.dirichlet.c*ones(ndf, 1);

% We consider given influx on the lower-left cell.

       injection_cells = 1; % Corresponds to lower-left cell for the Cartesian grid we have
                            % defined.
       bc.injection.cells = injection_cells;

     otherwise

       error('sim_case not recognized.')
   end
end

