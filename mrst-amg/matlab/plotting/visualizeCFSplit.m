function varargout = visualizeCFSplit(G, partitions, style, subset)
    if nargin == 2
        style = 'simple';
    end
    if nargin < 4
        subset = 1:numel(partitions);
    end
    
    level = nan(G.cells.num, 1);
    
    cellset = (1:G.cells.num)';
    for i = 1:numel(partitions)
        coarse = partitions{i}.coarse;
        
        level(cellset(coarse)) = i;
        cellset = cellset(coarse);
    end
    if nargout
        varargout{1} = level;
    else
        switch lower(style)
            case 'simple'
                plotCellData(G, level)
            case 'circle'
                wasHeld = ishold();
                hold on
                nl = numel(partitions);
                colors = (jet(nl));
                for i = subset
                    d = G.cells.centroids(level >= i, :);
                    X = d(:, 1);
                    Y = d(:, 2);
                    
                    if G.cells.num > 4096
                        mz = 4;
                    else
                        mz = 8;
                    end
                    if G.griddim == 2
                        plot3(X, Y, 0*X - 1,  'wo', 'MarkerFaceColor', colors(i, :), 'MarkerSize', mz);
                    else
                        Z = d(:, 3);
                        plot3(X, Y, Z, 'wo', 'MarkerFaceColor', colors(i, :), 'MarkerSize', mz);
                    end
                end
                
                if ~wasHeld
                    hold off
                end
            otherwise
                error('Unknown plotter');
        end
        
    end
end