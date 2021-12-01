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
B = 2*A
disp("A = ", A, "B = ", B)
C = A + B
D = A * B
disp("A + B = ", C, "A * B", D)
c = cond(A)
disp("conditionnement de A = ", c)

//EXO 7
eav = 0
ear = 0
ear2 = 0
for i = 1:1000
    A = rand(3,3)
    xex = rand(3,1)
    b = A * xex
    x_sol = A\b
    eav = eav + norm(xex - x_sol)/norm(xex)
    ear = ear + norm(b - A*x_sol)/(norm(A)*norm(x_sol))
    ear2 = ear2 + norm(b - A*x_sol)/norm(b)
end

disp("erreur avant moyen = ", eav/1000)
disp("erreur arrière (norme de A) moyen = ", ear/1000)
disp("erreur arrière (norme de b) myen = ", ear2/1000)

// for n = [100,1000,10000]
//     disp("n = ", n)
//     A = rand(n,n)
//     xex = rand(n,1)
//     b = A * xex
//     x_sol = A\b
//     err_av = norm(xex - x_sol)/norm(xex)
//     disp("erreur avant = ", err_av)
//     err_ar = norm(b - A*x_sol)/norm(b)
//     disp("erreur arrière = ", err_ar)
// end
