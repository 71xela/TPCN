function [X, nbiter, R] = jacobi_poisson1D(A,X,b,epsilon,nbiter)
    r = b-A*X
    R = [norm(r)]
    while norm(r) > epsilon
        X = X + 0.5*r
        r = b-A*X
        R = [R norm(r)]
        nbiter = nbiter + 1
    end
endfunction

n = 100
epsilon = 10**(-8)
T0 = -5.0
T1 = 5.0
v1 = -ones(1,n-1)
v2 = 2*ones(1,n)
grid = linspace(0,1,n+2)'
nbiter = 0

A = diag(v1, -1) + diag(v2, 0) + diag(v1, 1)
b = zeros(n,1)
b(1) = T0
b(n) = T1
EX_SOL = T0 + grid(2:n+1)*(T1 - T0)

X = zeros(n,1)

[X, nbiter, R] = jacobi_poisson1D(A,X,b,epsilon,nbiter)

//write("conv.txt", [R]);
disp(X)
disp(nbiter)


nbiter_tab = zeros(1,100)
for n=3:102
    disp("n = ", n)
    v1 = -ones(1,n-1)
    v2 = 2*ones(1,n)
    grid = linspace(0,1,n+2)'
    nbiter = 0

    A = diag(v1, -1) + diag(v2, 0) + diag(v1, 1)
    b = zeros(n,1)
    b(1) = T0
    b(n) = T1
    EX_SOL = T0 + grid(2:n+1)*(T1 - T0)

    X = zeros(n,1)

    [X, nbiter, R] = jacobi_poisson1D(A,X,b,epsilon,nbiter)
    nbiter_tab(n-2) = nbiter
end

// write("comp.txt", [nbiter_tab]);
