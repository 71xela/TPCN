exec("usolve.sci",-1)
exec("lsolve.sci",-1)

n = 10000
A = rand(n,n)
b = rand(n,1)
x = zeros(n,1)

L = tril(A)
U = triu(A)

tic()
[x] = usolve(U,b)
tU = toc()

tic()
[x] = lsolve(L,b)
tL = toc()

disp("Temps usolve : ", tU)
disp("Temps lsolve : ", tL)
