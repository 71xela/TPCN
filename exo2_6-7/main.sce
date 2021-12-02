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
minav = 0
maxav = 0
minar = 0
maxar = 0
minar2 = 0
maxar2 = 0
for i = 1:1000
    A = rand(3,3)
    xex = rand(3,1)
    b = A * xex
    x_sol = A\b
    av = norm(xex - x_sol)/norm(xex)
    ar = norm(b - A*x_sol)/(norm(A)*norm(x_sol))
    ar2 = norm(b - A*x_sol)/norm(b)

    if minav >= av then
        minav = av
    end

    if maxav <= av then
        maxav = av
    end

    if minar >= ar then
        minar = ar
    end

    if maxar <= ar then
        maxar = ar
    end

    if minar2 >= av then
        minar2 = ar2
    end

    if maxar2 <= av then
        maxar2 = ar2
    end

    eav = eav + av
    ear = ear + ar
    ear2 = ear2 + ar2
end

disp("erreur avant moyen = ", eav/1000)
disp("minav, maxav : ", minav, maxav)
disp("erreur arrière (norme de A) moyen = ", ear/1000)
disp("minar, maxar : ", minar, maxar)
disp("erreur arrière (norme de b) moyen = ", ear2/1000)
disp("minar2, maxar2 : ", minar2, maxar2)


for n = [100,1000] //[100,1000,10000]
    eav = 0
    ear = 0
    minav = 0
    maxav = 0
    minar = 0
    maxar = 0
    t = 0
    for i = 1:20
        A = rand(n,n)
        xex = rand(n,1)
        b = A * xex

        tic()
        x_sol = A\b
        t = t + toc()

        av = norm(xex - x_sol)/norm(xex)
        ar = norm(b - A*x_sol)/norm(b)

        if minav >= av then
            minav = av
        end

        if maxav <= av then
            maxav = av
        end

        if minar >= ar then
            minar = ar
        end

        if maxar <= ar then
            maxar = ar
        end

        eav = eav + av
        ear = ear + ar
    end
    disp("n = ", n)
    disp("erreur avant moyen = ", eav/20)
    disp("minav, maxav : ", minav, maxav)
    disp("erreur arrière (norme de b) moyen = ", ear/20)
    disp("minar, maxar : ", minar, maxar)
    disp("Temps fonction backslash : ", t/20)
end

//Une trentaine de secondes pour exécuter les lignes ci-après
n = 10000
disp("n = ", n)
A = rand(n,n)
xex = rand(n,1)
b = A * xex
tic()
x_sol = A\b
t = toc()
err_av = norm(xex - x_sol)/norm(xex)
disp("erreur avant = ", err_av)
err_ar = norm(b - A*x_sol)/norm(b)
disp("erreur arrière (norme de b) moyen = ", err_ar)
disp("Temps fonction backslash : ", t)
