function p = partitionAgglomerationLevels(A, coupleTol, checkDD, levelcount)
    [parts, Rs, Is] = deal(cell(levelcount, 1));
    
    for i = 1:levelcount
        % Partition at each level. Do diagonal dominance check at finest
        % level.
        parts{i} = partitionAgglomerationPairwise(A, coupleTol, checkDD*(i==1));
        [R, I] = unitInterpolator([], [], parts{i}, false);
        Rs{i} = R;
        Is{i} = I;
        A = R*A*I;
    end
    
    p = parts{end};
    for i = (levelcount-1):-1:1
        % Move up stack and store partition at each level
        I = Is{i};
        p = I*p;
    end
end