m = 100
n = 100
p = 100
A = rand(m,p)
B = rand(p,n)

getd(".")

tic()
[C] = matmat1b(A,B)
t1b= toc()

tic()
[C] = matmat2b(A,B)
t2b= toc()

tic()
[C] = matmat3b(A,B)
t3b= toc()

disp("Temps 1 boucle = ", t1b)
disp("Temps 2 boucle = ", t2b)
disp("Temps 3 boucle = ", t3b)
