function C = mtimes(A,B)
isAdic = strcmp(class(A),'dictionary');
isBdic = strcmp(class(B),'dictionary');
if isAdic && ~isBdic
    C = A.phi*B;
elseif isBdic && ~isAdic
    C = A*B.phi;
elseif isAdic && isBdic
    C = A.phi*B.phi;
end
end