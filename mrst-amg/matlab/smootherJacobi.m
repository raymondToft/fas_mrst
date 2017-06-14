function x = smootherJacobi(A, q, x, n, omega)
    D_inv = 1./diag(A);
    R = A - diag(diag(A));
    for i = 1:n
        x = omega*D_inv.*(q-R*x) + (1-omega)*x;
    end
end