function M = thresh_diagonal(M,n)


[n1,n2] = size(M);
M((n+1:n2) + (n:n2-1)*n1) = 0;



