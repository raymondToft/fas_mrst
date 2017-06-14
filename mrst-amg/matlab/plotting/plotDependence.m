function plotDependence(G, p)
    for i = 1:numel(p.dependence)
        dep = p.dependence{i};
        plotCellData(G, repmat(i, numel(dep), 1), dep, 'FaceAlpha', .2)
    end
    hold on
    gc = G.cells.centroids(p.C, :);
    X = gc(:, 1); Y = gc(:, 2);
    if G.griddim == 3
        Z = gc(:, 3);
        plot3(X, Y, Z, 'o', 'MarkerFaceColor', 'w');
    else
        plot(X, Y, 'o', 'MarkerFaceColor', 'w');
    end
end