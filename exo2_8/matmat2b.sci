function [C] = matmat2b(A,B)
    m = size(A)(1)
    n = size(B)(2)
    C = zeros(m,n)
    for i = 1:m
        for j = 1:n
            C(i,j) = A(i,:)*B(:,j) + C(i,j)
        end
    end
endfunction

