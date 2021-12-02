exec("exo3_2/usolve.sci",-1)
exec("exo3_2/lsolve.sci",-1)

for n = [10,20,50,70,100,200]
	disp("n = ", n)
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
        //xu = zeros(n,1)
        xl = zeros(n,1)
        L = tril(A)
		//U = triu(A)

		//bu = U*xex
		bl = L*xex

		//[xu] = usolve(U,bu)
        tic()
		[xl] = lsolve(L,bl)
        t = t + toc()

        av = norm(xex - xl)/norm(xex)
        ar = norm(bl - L*xl)/norm(bl)

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
    disp("erreur avant moyen = ", eav/20)
    disp("minav, maxav : ", minav, maxav)
    disp("erreur arrière (norme de b) moyen = ", ear/20)
    disp("minar, maxar : ", minar, maxar)
    disp("Temps descente : ", t/20)
end

n = 100
A = rand(n,n)
xex = rand(n,1)

U = triu(A)
b = U*xex

x = zeros(n,1)
[x] = usolve(U,b)
disp(abs(xex-x)./xex)
