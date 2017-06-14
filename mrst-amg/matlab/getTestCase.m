function [G, A, q, rock] = getTestCase(varargin)
    
    opt = struct('perm',    'homo',...
                 'gridaspect', 1, ...
                 'griddims', [], ...
                 'grid',    'cart', ...
                 'std',     10, ...
                 'speoffset', 0, ...
                 'stencil', 'tpfa', ...
                 'flow',    'bc', ...
                 'ordering', 'default', ...
                 'twistval',    0.03, ...
                 'layers',    1, ...
                 'dims',    5);
    opt = merge_options(opt, varargin{:});
    
    dims = opt.dims;
    if numel(dims) == 1
        dims = [1, 1]*2^dims;
    end
    
    if isempty(opt.griddims)
        opt.griddims = dims;
    end
    
    if ispc()
        dpath = 'D:\Jobb\data\';
    else
        dpath = '/data/';
    end
    
    switch lower(opt.grid)
        case 'cart'
            G = cartGrid(dims, opt.griddims.*opt.gridaspect);
        case 'twist'
            G = cartGrid(dims, opt.griddims.*opt.gridaspect);
            G = twister(G, opt.twistval);
        case 'l'
            G0 = cartGrid(dims, opt.griddims.*opt.gridaspect);
            ijk = gridLogicalIndices(G0);
            G = extractSubgrid(G0, ijk{1} < max(ijk{1})/2 | ijk{2} < max(ijk{2})/2);
            G.cartDims = G0.cartDims;
            clear G0
        case 'circle'
            G0 = cartGrid(dims, ones(size(dims)));
            G0 = computeGeometry(G0);
            
            X = G0.cells.centroids(:, 1);
            Y = G0.cells.centroids(:, 2);
            
            X = X - .5;
            Y = Y - .5;

            G = extractSubgrid(G0, sqrt(X.^2 + Y.^2) < .5);
            G.cartDims = G0.cartDims;
            clear G0
        case 'gullfaks'
            tmp = load(fullfile(dpath, 'gullfaks'));
            G = tmp.G;
        case 'norne'
            tmp = load(fullfile(dpath, 'norne'));
            G = tmp.G;
        otherwise
            error('Unknown grid')

    end
    
    mrstModule add gridtools coarsegrid
    switch lower(opt.ordering)
        case 'default'
        case 'cmk'
            G = reorderCMK(G, false);
        case {'red-black', 'rb'}
            G = reorderRedBlack(G);
        case 'amd'
            G = reorderAMD(G);
        case 'rcmk'
            G = reorderCMK(G, true);
    end
    
    try
    	G = mcomputeGeometry(G);
    catch
        G = computeGeometry(G);
    end
    
    
    bc = [];
    W = [];
    src = [];
    switch(lower(opt.perm))
        case {'homo', 'uniform'}
            rock.perm = ones(G.cells.num, 1);
            rock.poro = ones(G.cells.num, 1);
        case 'layers'
            if G.griddim == 3
                d = dims;
            else
                d = [1, dims];
            end
            rock.perm = logNormLayers(d, opt.layers, 'std', opt.std).*milli*darcy;
            
        case {'tarbert', 'ness'}
            mrstModule add spe10
            if strcmpi(opt.perm, 'tarbert')
                offset = 0 + opt.speoffset;
            else
                offset = 35 + opt.speoffset;
            end
            if numel(dims) == 2
                dims = [dims, 1];
            end
            rock = SPE10_rock(1:dims(1), 1:dims(2), (1:dims(3)) + offset);
            if G.griddim == 2
                rock.perm = rock.perm(:, 1);
            end
            rock.perm = convertFrom(rock.perm, milli*darcy);
            rock.poro(rock.poro < 0.01) = 0.01;
        case 'checker'
            require gridtools
            [ijk] = gridLogicalIndices(G);
            rock.perm = ones(G.cells.num, 1);
            for i = 1:G.griddim
                subs = ijk{i} < mean(ijk{i});
                rock.perm(subs) = rock.perm(subs).*(100 + rand(sum(subs), 1)*100);
            end
            rock.poro = rock.perm - min(rock.perm);
            rock.poro = rock.poro./max(rock.poro);
        case 'stripes'
            [i j] = gridLogicalIndices(G);
            rock.perm = ones(G.cells.num, 1);
            rock.poro = ones(G.cells.num, 1);
            
            rock.perm(mod(i, 2) == 0) = 0.1*rock.perm(mod(i, 2) == 0);
        case 'anisotropy'
            rock.perm = ones(G.cells.num, 2);
            rock.poro = ones(G.cells.num, 1);
            
            rock.perm(:, 2) = 0.1*rock.perm(:, 1);

        case 'jump2d'
            assert(G.griddim == 2);
            
            rock.perm = ones(G.cells.num, 2);
            rock.poro = ones(G.cells.num, 1);
            [i j] = gridLogicalIndices(G);
            
            X = normalize(G.cells.centroids(:, 1));
            Y = normalize(G.cells.centroids(:, 2));
            
            rock.perm(X >= .65 & X <= .95 & Y >= .05 & Y <= .65, 2) = 100;
            rock.perm(X >= .25 & X <= .45 & Y >= .25 & Y <= .45, 1) = 100;
            rock.perm(X >= .05 & X <= .25 & Y >= .65 & Y <= .95, :) = 100;
            
            
            if strcmpi(opt.flow, 'none')
                srccells = find(X >= .05 & X <= .25 & Y >= .65 & Y <= .95);
                src = addSource(src, srccells, ones(numel(srccells), 1));
                bc = pside(bc, G, 'ymax', 0);
            end
% 
%             
%             mi = G.cartDims(1);
%             mj = G.cartDims(2);
%             rock.perm(i > mi/4 & i < 2*mi/4 & j > mj/4 & j < 2*mj/4, 2) = 1e-3;
%             rock.perm(i > 2*mi/3 & i < 15*mi/16 &...
%                       j > mj/16   & j < 2*mj/3, 2) = 1e-3;
% 
%             
        case 'grid'
            rock = tmp.rock;
        otherwise
            error('Unknown perm case');
    end
    
    fluid = initSingleFluid('rho', 1, 'mu', 1);
    sol = initResSol(G, 0);

    switch lower(opt.flow)
        case 'bc'
            bc = pside(bc, G, 'xmin', 0);   
            bc = pside(bc, G, 'xmax', 1);
        case 'none'
            
        case 'wells'
            W = verticalWell(W, G, rock, 1, 1, [], 'val', 1);
            W = verticalWell(W, G, rock, dims(1), dims(2), [], 'val', 0);
        case 'src'
            src = addSource(src, 1, 1);
            src = addSource(src, G.cells.num, -1);
        otherwise
            error()
    end
    
    gravity off
    
    switch lower(opt.stencil)
        case 'mpfa'
            mrstModule add mpfa
            T = computeMultiPointTrans(G, rock);
            solver = @(varargin) incompMPFA(varargin{:});
        case 'tpfa'
            T = computeTrans(G, rock);
            solver = @(varargin) incompTPFA(varargin{:});
        otherwise
            error('Unknown solver');
    end
    sol = solver(sol, G, T, fluid, 'bc', bc, 'wells', W, 'src', src, 'matrixOutput', true);
    
    subs = 1:G.cells.num;
    A = sol.A(subs, subs);
    q = sol.rhs(subs);
    if isempty(bc) && isempty(W) && isempty(src)
        A(1,1) = A(1,1)/2;
    end
end

function d = normalize(d)
    d = d - min(d);
    d = d./max(d);
end
