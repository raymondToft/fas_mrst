%%
% 
% perlOneLineDescription(Setup the system)
%  
%  At the beginning of the simulation, we compute the non-dynamic variables such as
%
% * the pore volumes, |pv|,
% * the transmissibilities, |T|, 
% * the discrete differential operators, |div| and |grad|,
% * the gravity force contribution, |dz|,
%
%%

function s = setupSystem(G, rock, bc, param)

   cf = G.cells.faces(:,1);
   nf = G.faces.num;
   nc = G.cells.num;
   
   %% 
   % Compute pore volumes
   % 
   
   s.pv = poreVolume(G, rock);
   s.poro = rock.poro;

   %%
   % Compute the half, and then the full, transmissibilities.
   %
   
   T = computeTrans(G, rock);
   cf = G.cells.faces(:,1);
   nf = G.faces.num;
   T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
   s.T = T;

   %% 
   % Set up the discrete divergence operator, |div|. It sums up all signed faces'
   % values in each cell.
   % 
   
   N = double(G.faces.neighbors);
   index = 1:nf';
   faces1 = N(:, 1) ~= 0;
   faces2 = N(:, 2) ~= 0;
   C1  = sparse(index(faces1), N(faces1, 1), ones(nnz(faces1),1), nf, nc);
   C2  = sparse(index(faces2), N(faces2, 2), ones(nnz(faces2),1), nf, nc);
   C = C1 - C2;
   s.div  = @(x)C'*x;

   %%
   % Set up the discrete gradient operator, |grad|. It is an operator from cell values to
   % face values. We compute the differences of cell values across each face. 
   
   index = 1:nf;
   interior = prod(N, 2)~=0;

   C1_interior = sparse(index(interior), N(interior, 1), ones(nnz(interior), 1), nf, nc);
   C2_interior = sparse(index(interior), N(interior, 2), ones(nnz(interior), 1), nf, nc);

   
   %%
   % Compute the boundary contribution to the gradient operator. They corresponds to the
   % external faces where Dirichlet conditions are given. We are careful to use the
   % correct signs.
   % 

   is_dirichlet_faces1 = N(bc.dirichlet.faces, 1) ~= 0;
   is_dirichlet_faces2 = ~is_dirichlet_faces1;

   dirichlet_faces1 = bc.dirichlet.faces(is_dirichlet_faces1);
   dirichlet_faces2 = bc.dirichlet.faces(is_dirichlet_faces2);

   C1_exterior = sparse(index(dirichlet_faces1), ...
                        N(dirichlet_faces1, 1), ...
                        ones(numel(dirichlet_faces1), 1), nf, nc);
   C2_exterior = sparse(index(dirichlet_faces2), ...
                        N(dirichlet_faces2, 2), ...
                        ones(numel(dirichlet_faces2), 1), nf, nc);

   %%
   % The gradient operator is the sum of internal and boundary contributions. For
   % boundary faces, if they are Dirichlet faces, then we use the Dirichlet value to
   % compute the difference, otherwise the gradient is set to zero, which corresponds
   % to a no flow condition.
   %
   
   C = C1_interior + C1_exterior - (C2_interior + C2_exterior);

   pressure_bc = sparse(nf, 1);
   pressure_bc(dirichlet_faces1) = - bc.dirichlet.pressure(is_dirichlet_faces1);
   pressure_bc(dirichlet_faces2) = + bc.dirichlet.pressure(is_dirichlet_faces2);

   s.p_grad = @(p)(C*p + pressure_bc);  

   s.grad = @(val, bc_val)(grad(val, bc_val, nf, C, ...
                                is_dirichlet_faces1, dirichlet_faces1, ... 
                                is_dirichlet_faces2, dirichlet_faces2));

   %%
   % Set up the gravity term.
   %
   
   z = G.cells.centroids(:, 3);
   fz = G.faces.centroids(:, 3);
   s.dz = s.grad(z, fz);

   %% 
   % We define the function that return values on faces using cells' values and upwind
   % directions. For boundary faces when the upwind direction require to take value
   % from the outside, if the face is a Dirichlet face, then we take the Dirichlet
   % value otherwise we take the value of the neighboring cell.
   %
   
   s.upwindFaceValues = @(flag, val, bc_val) upwindFaceValues(flag, val, bc_val, N, interior, ...
                                                     dirichlet_faces2, dirichlet_faces1,  bc, nf, ...
                                                     nc);

   %% 
   % We define the function that return values on faces defined as the average of
   % neighboring cells. For boundary faces, we use the value of the neighboring cell if
   % it is not a Dirichlet cell or the average of this value and the Dirichlet face
   % value. The Dirichlet face values are specified with |bc_bal|.
   %
   
   s.averageFaceValues = @(val, bc_val) averageFaceValues(val, bc_val, N, interior, ...
                                                     dirichlet_faces2, dirichlet_faces1,  bc, nf, ...
                                                     nc);
   
   s.N = N;
   s.G = G;

end


%%
% <html><h2>Helper functions</h2></html>
%


%%
% perlInsertHtml("<h3>Helper function for the |grad| operator</h3>")
%

function dval = grad(val, bc_val, nf, C, ...
                     is_dirichlet_faces1, dirichlet_faces1, ...
                     is_dirichlet_faces2, dirichlet_faces2)

   signed_bc_val = sparse(nf, 1);
   signed_bc_val(dirichlet_faces1) = - bc_val(is_dirichlet_faces1);
   signed_bc_val(dirichlet_faces2) = + bc_val(is_dirichlet_faces2);
   dval = C*val + signed_bc_val;

end


%%
% perlInsertHtml("<h3>Helper function for the |upwindFaceValues| operator</h3>")
%


function up_val = upwindFaceValues(flag, val, bc_val, N, interior, dirichlet_faces2, dirichlet_faces1, ...
                                     bc, nf, nc)
   index        = (1:nf)';
   upCell       = N(:, 2);
   upCell(flag) = N(flag, 1);

   %%
   % On the interior cell we use upwind
   %
   
   Mint = sparse(index(interior), upCell(interior), 1, nf, nc);

   logical_dirichlet_faces1 = zeros(nf, 1);
   logical_dirichlet_faces1(dirichlet_faces1) = 1;
   logical_dirichlet_faces1 = logical(logical_dirichlet_faces1);
   logical_dirichlet_faces2 = zeros(nf, 1);
   logical_dirichlet_faces2(dirichlet_faces2) = 1;
   logical_dirichlet_faces2 = logical(logical_dirichlet_faces2);

   external_faces1 = N(:,2)==0;
   external_faces2 = N(:,1)==0;

   %%
   % On the exterior faces where no Dirichlet conditions are given we take the value given in
   % the interior cell.
   %
   
   Mext1 = sparse(index(external_faces1 & ~logical_dirichlet_faces1), ...
                  N(external_faces1 & ~ logical_dirichlet_faces1, 1), 1, nf, nc);
   Mext2 = sparse(index(external_faces2 & ~logical_dirichlet_faces2), ...
                  N(external_faces2 & ~ logical_dirichlet_faces2, 2), 1, nf, nc);

   %%
   % On the Dirichlet boundary cells we use upwind, taking the values from boundary
   % conditions when needed We assume flag is logical
   %
   
   assert(islogical(flag), 'Upwind indices must be given as logical');
   Mdir1 = sparse(index(flag & logical_dirichlet_faces1), ...
                  N((flag & logical_dirichlet_faces1), 1), 1, nf, nc);
   Mdir2 = sparse(index(~flag & logical_dirichlet_faces2), ...
                  N((~flag & logical_dirichlet_faces2), 2), 1, nf, nc);

   M = Mint + Mext1 + Mext2 + Mdir1 + Mdir2;

   %%
   % Values of saturation from Dirichlet boundary conditions.
   %
   
   dval_all = sparse(nf, 1);
   dval_all(bc.dirichlet.faces) = bc_val; 
   dval = sparse(nf, 1);
   dval(flag & logical_dirichlet_faces2)  = dval_all( flag & logical_dirichlet_faces2);
   dval(~flag & logical_dirichlet_faces1) = dval_all(~flag & logical_dirichlet_faces1);
   
   up_val = M*val + dval;
   
end


%%
% perlInsertHtml("<h3>helper function for the |averageFaceValues| operator</h3>")
% 
%


function av_val = averageFaceValues(val, bc_val, N, interior, dirichlet_faces2, dirichlet_faces1, ...
                                    bc, nf, nc)
   index        = (1:nf)';

   % On the interior cell we use average

   Mint1 = sparse(index(interior), N(interior, 1), 0.5, nf, nc);
   Mint2 = sparse(index(interior), N(interior, 1), 0.5, nf, nc);

   logical_dirichlet_faces1 = zeros(nf, 1);
   logical_dirichlet_faces1(dirichlet_faces1) = 1;
   logical_dirichlet_faces1 = logical(logical_dirichlet_faces1);
   logical_dirichlet_faces2 = zeros(nf, 1);
   logical_dirichlet_faces2(dirichlet_faces2) = 1;
   logical_dirichlet_faces2 = logical(logical_dirichlet_faces2);

   % For cells at a boundary which is not dirichlet, we use the cell value
   exterior = ~(logical_dirichlet_faces1 | logical_dirichlet_faces2 | interior);
   Mext = sparse(index(exterior), N(exterior, 1) + N(exterior, 2), 1, nf, nc);

   % On the Dirichlet boundary cells we take average of cell value and dirichlet boundary value

   Mdir1 = sparse(index(logical_dirichlet_faces1), N(logical_dirichlet_faces1, 1), 0.5, nf, nc);
   Mdir2 = sparse(index(logical_dirichlet_faces2), N(logical_dirichlet_faces2, 2), 0.5, nf, nc);

   M = Mint1 + Mint2 + Mext + Mdir1 + Mdir2;

   % Values of saturation from Dirichlet boundary conditions.
   dval = sparse(nf, 1);
   dval(bc.dirichlet.faces) = 0.5*bc_val; 
   
   av_val = M*val + dval;
   
end