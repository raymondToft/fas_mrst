function plotCoarseData(coarsegrids, index, data, varargin)
    G_fine = coarsegrids{1}.parent;
    
    if index == 0
        dataGrid = G_fine;
    else
        dataGrid = coarsegrids{index};
    end
    assert(any(size(data) == 1));
    if mod(numel(varargin), 2)
        subs = varargin{1};
        if islogical(subs)
            subs = find(subs);
        end
        varargin = varargin(2:end);
        
        tmp = nan(dataGrid.cells.num, 1);
        tmp(subs) = data;
        data = tmp;
        clear tmp
    end
    
    if index
        for i = index:-1:1
            data = data(coarsegrids{i}.partition);
        end
    end
    plotCellData(G_fine, data, varargin{:});
end