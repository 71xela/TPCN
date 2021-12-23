function [A] = ldltsym(A)
    n = (size(A))(1)
    v = zeros(n,1)

    for j = 1:n
        for i = 1:j-1
            v(i) = A(j,i)*A(i,i)
        end
        A(j,j) = A(j,j) - A(j,1:j)*v(1:j)
        A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n,1:j)*v(1:j))/A(j,j)
    end
endfunction

function [A] = lu1b(A)
    n = size(A)(1)
    for k = 1:n-1
        A(k+1:n,k) = A(k+1:n,k)/A(k,k)
        A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - A(k+1:n,k)*A(k,k+1:n)
    end
endfunction


//main
A = [[10,20,30];[20,45,80];[30,80,171]]

[A] = ldltsym(A)
disp(A)

t = 0
for n=[10,20,50,70,100,150,200,500,700,1000]
    disp("n = ", n)
    A = rand(n,n) + eye(n,n)*n

    for i=1:20
        tic()
        [A] = ldltsym(A)
        t = t + toc()
    end
    disp("Temps factorisation LDL^T : ", t/20)
end

t = 0
for n=[10,20,50,70,100,150,200,500,700,1000]
    disp("n = ", n)
    A = rand(n,n) + eye(n,n)*n

    for i=1:20
        tic()
        [A] = lu1b(A)
        t = t + toc()
    end
    disp("Temps factorisation LU : ", t/20)
end