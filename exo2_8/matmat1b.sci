function [C] = matmat1b(A,B)
    m = size(A)(1)
    n = size(B)(2)
    C = zeros(m,n)
    for i = 1:m
        C(i,:) = A(i,:)*B + C(i,:)
    end
endfunction

