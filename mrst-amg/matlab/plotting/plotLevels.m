function plotLevels(levels)
    if isfield(levels, 'level')
    	levels = levels.level;
    end
    hold on
    plot(max(levels) - levels, '-o', 'linewidth', 2, 'markerfacecolor', 'r')
    
    exact = find(levels == max(levels));
    plot(exact, zeros(size(exact)), 'o', 'markerfacecolor', 'w');
    
    levels = unique(levels);
    set(gca, 'YTickMode', 'manual')
    set(gca, 'YTick', levels);
    levelnames = arrayfun(@num2str, levels, 'unif', false);
    
    levelnames{1} = '(Coarsest)';
    levelnames{end} = '(Finest)';
    set(gca, 'Yticklabel', levelnames)
    legend('Smoother', 'Exact solver')
    grid on
    axis tight
%     title('Cycle structure')
end