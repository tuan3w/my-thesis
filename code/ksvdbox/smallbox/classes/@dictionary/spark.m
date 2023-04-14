function varargout = Spark(obj)
% Calculates the minimum number of linearly dependent atoms in the matrix A
% WARNING: this function computes a combinatorial function, use only if the
% size of the problem is small (i.e. <20)
if nargout <= 1
    A = obj.phi;
    k = size(A,2);
    if k>20
        warning('This function computes a combinatorial function, use only if thesize of the problem is small (i.e. <20).');
        choice = input('The calculation of spark will take a long time... do you wish to continue anyway (y/n)','s');
        if strcmpi( choice,'n')
            return
        end
    end
    sigma = 2;
    while true
        P = nchoosek(1:k,sigma);
        for i=1:size(P,1)
            r = rank(A(:,P(i,:)));
            if r<sigma
                varargout{1} = sigma;
                return
            end
        end
        sigma = sigma + 1;
        if sigma==k
            varargout{1} = inf;
            return
        end
    end
else
    %% TODO: calculate lower and upper bounds on the spark
    varargout{1} = 2;
    varargout{2} = inf;
end
