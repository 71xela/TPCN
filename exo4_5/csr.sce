function [AX] = csrsolve(AX, AA, JA, IA, V, n)
    for i=1:n
        AX(i) = AA(IA(i):IA(i+1)-1) * V(JA(IA(i):IA(i+1)-1))
    end
endfunction

function [B] = matvec1b(B, A, V, n)
    for i=1:n
        for j=1,n
            B(i) = B(i) + A(i,j) * V(j)
        end
    end
endfunction

A = [[15, 0, 0, 22, 0, -15, 0, 0];
    [0, 11, 3, 0, 0, 0, 2, 0];
    [0, 0, 0, -6, 0, 0, 0, 0];
    [0, 0, 0, 0, 0, 0, 0, 0];
    [91, 0, 0, 0, 0, 0, 25, 7];
    [0, 0, 28, 0, 0, 0, 0, -2]]

n = 6
AX = zeros(n,1)
B = zeros(n,1)

// Stockage CSR
AA = [15, 22, -15, 11, 3, 2, -6, 91, 25, 7, 28, -2]
JA = [1, 4, 6 ,2 ,3, 7, 4, 1, 7, 8, 3, 8]
IA = [1, 4, 7, 8, 8, 11, 13]

// Vecteur
V1 = ones(8,1)
V2 = [1, 0, 0, 1, 0, 0, 1, 0]'

tic()
[B] = matvec1b(B, A, V1, n)
t = toc()
disp("Matrice-v1 BLAS : ", B)
disp("Temps BLAS : ", t)
tic()
[AX] = csrsolve(AX, AA, JA, IA, V1, n)
t = toc()
disp("Matrice-v1 CSR : ", AX)
disp("Temps CSR : ", t)

AX = zeros(n,1)
B = zeros(n,1)

tic()
[B] = matvec1b(B, A, V2, n)
t = toc()
disp("Matrice-v2 BLAS : ", B)
disp("Temps BLAS : ", t)
tic()
[AX] = csrsolve(AX, AA, JA, IA, V2, n)
t = toc()
disp("Matrice-v2 CSR : ", AX)
disp("Temps CSR : ", t)
