exec("mylu3b.sci", -1)
exec("mylu1b.sci", -1)

n = 100
A = rand(n,n)

tic()
[L,U] = mylu3b(A)
t = toc()
disp("Temps 3 boucles = ", t)
disp(sum(A-L*U))

tic()
[L,U] = mylu1b(A)
t = toc()
disp("Temps 1 boucles = ", t)
disp(sum(A-L*U))
