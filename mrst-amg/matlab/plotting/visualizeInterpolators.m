function d = visualizeInterpolators(G, partitions)
    nl = numel(partitions);
    cellset = (1:G.cells.num)';
    
    d = struct();
    part_old = (1:G.cells.num)';
    
    for i = 1:nl
        I = partitions{i}.I;
        [~, part] =  max(I,[],2);
        
        name = ['level_', num2str(i)];
        
        d.(name).I = I(part_old, :);
        d.(name).partition = part_old;

        part_old = part(part_old);
        
        coarse = partitions{i}.coarse;
        cellset = cellset(coarse);

    end
end