[nx,ny,nz] = deal( 24,  24, 2);
[Dx,Dy,Dz] = deal(500, 500, 50);
grid = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
grid = computeGeometry(grid);

perm_range = [0.2 0.3];
gauss_filter_size = [3 3 3];
std = 2.5;

p = gaussianField(grid.cartDims, perm_range, gauss_filter_size, std);
K = p.^3. * (1e-5)^2./(0.81 * 72 * (1-p).^2);
rock = makeRock(grid, K(:), p (:));

[min(K), max(K)]


figure;
plotCellData (grid , convertTo ( rock.perm , milli * darcy ));
colorbar (); axis equal tight ; view (3);
  