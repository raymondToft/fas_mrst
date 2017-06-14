function x = smootherGS(A, q, x, n, up)
    if nargin == 4
        up = false;
    end
    
    if up
        % U then L
        L = tril(A, 0);
        U = triu(A, 1);
    else
        % L then U
        L = triu(A, 0);
        U = tril(A, -1);
    end
    for i = 1:n
        x = L\(q - U*x);
    end
end