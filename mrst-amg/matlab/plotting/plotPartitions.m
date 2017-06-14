function s = plotPartitions(G, partitions, simple)
    if nargin == 2
        simple = false;
    end
    require mrst-gui
    np = numel(partitions);
    colors = lines(np);
    s = struct('dummy', []);
    parts = nan(G.cells.num, np);
    
    doPlot = ~nargout;
    for i = 1:np
        part = partitions{i}.partition;
        for j = (i-1):-1:1
            p = partitions{j};
            [R, I] = unitInterpolator([], [], p, false);
            part = I*part;
        end
        if simple && doPlot
            c = colors(i, :);
            if G.griddim == 2
                outlineCoarseGrid(G, part, c, 'linestyle', '--', 'linewidth', 2);
            else
                outlineCoarseGrid(G, part, 'FaceColor', c,...
                                           'EdgeColor', min(1.1*c, 1),...
                                           'FaceAlpha', .1,...
                                           'EdgeAlpha', .2);
            end
        else
            s.(['Partition_', num2str(i)]) = mod(part, 13) + part/max(part);
            parts(:, i) = part;
        end
    end
    s.parts = parts;
    if ~simple && doPlot
        plotToolbar(G, s, 'EdgeColor', 'w', 'EdgeAlpha', .1);
    end
end