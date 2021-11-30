function [ out ] = FECTP( M, p, numIter)
    ## This functions returns the length(M) + 1-sized array [eigenvector, eigenvalue] that solves the eigenproblem of a hermitian  matrix M
    ## p serves as a guess for the eigenvalue.
    ## Based on this very cool algorithm https://en.wikipedia.org/wiki/Rayleigh_quotient_iteration
    n = length(M);
    w = rand(n, 1);
    u = rand(n, 1);
    w = u;
    u(1) = 0;
    u(end) = 0;
    PAJ = -p * speye(n);
    for s = 1 : numIter
        u = (M + PAJ)\u;
        u = u / norm(u);
        %% Optional code to include different boundary conditions, experimental.
        %u(1) = 0;
        %u(end) = 0;
    end
    out = [u', p + sum(w.*u)/sum(w.*((M + PAJ)\u))];
end
