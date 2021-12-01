//EXO 6
x = 1:4
y = (4:-1:1)'
disp("x = ", x, "y = ", y)
z = x + y'
s = x * y
disp("z = ", z, "s = ", s)
nx = size(x)
ny = size(y)
disp("taille x = ", nx, "taille y = ", ny)
norme2x = norm(x)
disp("||x|| = ", norme2x)
A = [1:4; 5:8; 9:12]
A_T = A'
disp("A = ", A, "Transposée de A = ", A_T)
A = [1:3; 4:6; 7:9]
B = A'
disp("A = ", A, "B = ", B)
C = A + B
D = A * B
disp("A + B = ", C, "A * B", D)
c = cond(A)
disp("conditionnement de A = ", c)

//EXO 7
A = rand(3,3)
disp("A = ", A)
xex = rand(3,1)
disp("xex = ", xex)
b = A * xex
disp("b = ", b)
x_sol = A\b
disp("Solution de Ax = b -->", x_sol)
err_av = norm(xex - x_sol)/norm(xex)
disp("erreur avant = ", err_av)
err_ar = norm(b - A*x_sol)/(norm(A)*norm(x_sol))
disp("erreur arrière = ", err_ar)

for n = [100,1000,10000]
    disp("n = ", n)
    A = rand(n,n)
    xex = rand(n,1)
    b = A * xex
    x_sol = A\b
    err_av = norm(xex - x_sol)/norm(xex)
    disp("erreur avant = ", err_av)
    err_ar = norm(b - A*x_sol)/(norm(b))
    disp("erreur arrière = ", err_ar)
end
