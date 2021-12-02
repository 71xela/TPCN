exec("exo3_3/gausskij3b.sci",-1)

for n=[10,20,50,70,100]
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
	    b = A * xex
	    x_sol = zeros(n,1)

		tic()
		[x_sol] = gausskij3b(A,b)
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
	disp("erreur avant moyen = ", eav/20)
	disp("minav, maxav : ", minav, maxav)
	disp("erreur arrière (norme de b) moyen = ", ear/20)
	disp("minar, maxar : ", minar, maxar)
	disp("Temps Gauss+usolve : ", t/20)
end
