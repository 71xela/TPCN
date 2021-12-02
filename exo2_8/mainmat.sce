exec("exo2_8/matmat1b.sci",-1)
exec("exo2_8/matmat2b.sci",-1)
exec("exo2_8/matmat3b.sci",-1)

for n=[10,20,50,70,100,150,200]
	disp("n = ", n)
	A = rand(n,n)
	B = rand(n,n)
	t1b = 0
	t2b = 0
	t3b = 0

	for i=1:10
		tic()
		[C] = matmat1b(A,B)
		t1b = t1b + toc()
	end

	for i=1:10
		tic()
		[C] = matmat2b(A,B)
		t2b = t2b + toc()
	end

	for i=1:5
		tic()
		[C] = matmat3b(A,B)
		t3b = t3b + toc()
	end

	disp("Temps 1 boucle = ", t1b/10)
	disp("Temps 2 boucle = ", t2b/10)
	disp("Temps 3 boucle = ", t3b/5)
end
