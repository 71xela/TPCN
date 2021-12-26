function [X, nbiter, R] = richardson_poisson1D(A,X,b,epsilon,nbiter, alpha)
    r = b-A*X
    R = [norm(r)]
    while norm(r) > epsilon
        X = X + alpha*r
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

//vp = spec(A)
vp = zeros(1,n)
for i=1:n
    vp(i) = 2-2*cos(%pi*i/(n+1))
end
lmin = vp(1)
lmax = vp(n)

alpha = 2/(lmax+lmin)
disp("alpha = ", alpha)

X = zeros(n,1)

[X, nbiter, R] = richardson_poisson1D(A,X,b,epsilon,nbiter,alpha)

// mod = string("w")
// fd = mopen("conv49.txt", mod)
// for i=R
//     mfprintf(fd, "%e\n",i)
// end
// mclose(fd)

disp(X)
disp(nbiter)
