function partition = standardCoarsening(A, coupleTol, varargin)

    assert((coupleTol <= 1 && coupleTol >= 0) || isinf(coupleTol))
    % Transform to positive, remove diagonal
    A_offdiag = abs(A - diag(diag(A)));
    connStr = full(max(A_offdiag, [], 2));
    
    % Undecided variables
    U = true(size(A,1),1);
    
    % Fine / coarse split
    [F, C] = deal(~U);
    
    % Get row and column indices for all connections
    [rowi, coli, conn] = find(A_offdiag);
    
    if numel(varargin)
        % Strict tolerance
        isStrong0 = conn >= varargin{1};
    else
        isStrong0 = conn > coupleTol*connStr(rowi);
    end
    isStrong = isStrong0;
    isStrongFine = zeros(size(isStrong));
    
    measure =  accumarray(rowi, isStrong);
    iter = 0;
    while any(measure > 0)
        % The next coarse variable is the globally largest value
        [v, nextCoarse] = max(measure); %#ok Backwards compatability
        
        % The next fine variables are those connected to the current coarse
        nextFine = rowi(coli == nextCoarse & isStrong);
        
        % Flag the updated values
        changed = false(size(F));
        changed(vertcat(nextFine, nextCoarse)) = true;
        update = changed(coli);
        % Subset corresponding to possibly updated values - those where the
        % column index is a pointer to a updated node
        subs = rowi(update);
        
        % Calculate the old measure for the subset
        oldMeasure =  accumarray(subs, isStrong(update)) +...
                    2*accumarray(subs, isStrongFine(update));
        
        % Update the maps
        F(nextFine) = true;
        C(nextCoarse) = true;
        isStrong(nextCoarse) = 0;
        isStrongFine(nextFine) = isStrong(nextFine);
        
        % Compute the new measured value for the update subset
        newMeasure =  accumarray(subs, isStrong(update)) + ...
                    2*accumarray(subs, isStrongFine(update));
        
        target = 1:max(subs);

        
        % Update the measure
        measure(target) = measure(target) - oldMeasure + newMeasure;
        measure(F | C) = 0;

        if 0
            figure(1);
            clf
            plotCellData(evalin('base', 'G'), measure)
            plotGrid(evalin('base', 'G'), F, 'FaceColor', 'none', 'EdgeColor', 'r', 'LineWidth', 2)
            plotGrid(evalin('base', 'G'), C, 'FaceColor', 'g', 'EdgeColor', 'g', 'LineWidth', 2)
            colorbar
            if isnan(c)
                c = [min(measure), max(measure)];
            else
                caxis(c)
            end
            unique(subs)
            pause();
            nextCoarse
        end
        iter = iter + 1;
        if mod(iter, 100) == 0
            dispif(mrstVerbose, '%d of %d\n', iter, size(A, 1));
        end
    end
    F = ~C;
   
    coarse = find(C);
    nc = numel(coarse);
    coarseConn = cell(nc, 1);
    for i = 1:nc
        coarseConn{i} = [coarse(i); unique(coli(rowi == coarse(i)))];
    end
    partition = struct('coarse',        C, ...
                       'fine',          F & ~C, ...
                       'strongConn',    isStrong0);
    partition.strongCoarseConnected = coarseConn;
end

%         % Sum of strong connections *TO* cell i
%         S =  accumarray(coli, isStrong);
%         % Sum of strong connections *FROM* cell i
%         St = accumarray(rowi, isStrong);
