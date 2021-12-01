exec("gausskij3b.sci",-1)

n=100
A = rand(n,n)
b = rand(n,1)
x = zeros(n,1)

tic()
[x] = gausskij3b(A,b)
t = toc()

disp("Temps Gauss+usolve : ", t)
