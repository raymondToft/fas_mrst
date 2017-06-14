function [x, meta] = multigridCycle(A, q, x, v1, v2, ...
                                    cycleindex, partitions,...
                                    smoother, getOperators, ...
                                    varargin)
    if numel(varargin) == 0
        % Instansiate meta and assume level to be 0
        meta = struct('level',  [], ...
                      'defect', []);
        level = 0;
    else
        meta = varargin{1};
        level = varargin{2};
    end
    
    if isempty(partitions)
        x = A\q;
        meta = storeLevel(meta, level, A, q, x);
        return
    end
    
    if v1
        % presmoother
        x = smoother(A, q, x, v1);
        meta = storeLevel(meta, level, A, q, x);
    end
    
    % Calculate defect and send to the coarser grid
    defect = q - A*x;
    
    [R, I] = getOperators(A, q, partitions{1});
    defect_coarse = R*defect;
    if isfield(partitions{1}, 'A')
        A_coarse = partitions{1}.A;
    else
        A_coarse = R*A*I;
    end
    correct_coarse = zeros(size(defect_coarse));
    
    iter = cycleindex;
    if numel(iter) == 1
        iter = cycleindex*ones(1, iter);
    end
    for i = iter
        [correct_coarse, meta] = multigridCycle(A_coarse, defect_coarse, correct_coarse,...
                                        v1, v2, ...
                                        i, partitions(2:end), ...
                                        smoother, getOperators,...
                                        meta, level + 1);
        if numel(partitions) == 1;
            % We are at lowest level, the correction is exact
            break
        end
    end
    x = x + I*correct_coarse;
    
    if v2
        % postsmoother
        x = smoother(A, q, x, v2);
        meta = storeLevel(meta, level, A, q, x);
    end
end

function meta = storeLevel(meta, level, A, q, x)
    def_norm = norm(A*x - q, inf)/norm(q, inf);
    meta.level = [meta.level; level];
    meta.defect = [meta.defect; def_norm];
end