function [L,U,P] = mylu(A)
    n = size(A)(1)
    for k = 1:n-1
        [piv,ind] = max(abs(A(k:n,k)))
        ind = k-1+ind
        q = row(1,ind)
        if ind ~= k then
            new = A(ind,:)
            A(ind,:) = A(k,:)
            row(1,ind) = row(1,k)
            row(1,k) = q
            A(k,:) = new
        end
        A(k+1:n,k) = A(k+1:n,k)/A(k,k)
        A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - A(k+1:n,k)*A(k,k+1:n)
    end
    L = tril(A)
    for i = 1:n
        L(i,i) = 1
    end
    U = triu(A)
endfunction
