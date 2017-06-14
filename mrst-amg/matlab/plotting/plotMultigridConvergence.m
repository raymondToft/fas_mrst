function plotMultigridConvergence(meta)
    meta = meta(~cellfun(@isempty, meta));
    set(gcf, 'renderer', 'painters')
    subplot(2, 2, 1)
    plotLevels(meta{1});
    subplot(2, 2, 2)
    plotFineConvergence(meta)
    
    subplot(2, 2, 3:4)
    plotSubconvergence(meta)
   
    
end


function plotFineConvergence(meta)
    res = cellfun(@(x) x.defect(end), meta);

    semilogy(res, 'linewidth', 2);
    axis tight
    title('Fine scale convergence')
end

function plotSubconvergence(meta)
    levels = cellfun(@(x) x.level, meta, 'unif', false);
    levels = vertcat(levels{:});
    
    defect = cellfun(@(x) x.defect, meta, 'unif', false);
    defect = vertcat(defect{:});
    
    l = unique(levels);
    nl = numel(l);
    nd = numel(defect);
    
    data = nan(nd, nl);
    
    for i = 1:numel(l)
        subs = levels == l(i);
        data(subs, i) = defect(subs);
    end
    
    plot(data, '.');
    axis tight
    grid on
    legend(arrayfun(@(x) ['Level: ' num2str(x)], l, 'unif', false), ...
                                        'location', 'eastoutside');
    title('Convergence, levels')
end