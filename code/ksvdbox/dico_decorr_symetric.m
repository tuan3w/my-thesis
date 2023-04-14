function dico = dico_decorr_symetric(dico, mu, B, centroids, nclass)
    %DICO_DECORR decorrelate a dictionary
    %   Parameters:
    %   dico: the dictionary
    %   mu: the coherence threshold
    %   amp: the amplitude coefficients, only used to decide which atom to
    %   project
    %
    %   Result:
    %   dico: a dictionary close to the input one with coherence mu.
    
    eps = 1e-6; % define tolerance for normalisation term alpha
    
    % convert mu to the to the mean direction
    theta = acos(mu  - 0.1)/2;
    ctheta = cos(theta);
    stheta = sin(theta);
    PMat = cell(nclass, nclass);
    for i=1:nclass
        for j = 1:nclass
            %D_i = normcols([centroids{i}, centroids{j}]);
            p1 = centroids{i} + centroids{j};
            p2 = centroids{i} - centroids{j};
            D_i = [p1, p2];
            %D_i = D_i ;
            PMat{i}{j} = D_i * inv(D_i'*D_i + 0.0000001* eye(size(D_i,2))) *D_i';
            PMat{j}{i} = PMat{i}{j};
        end
    end
    %pause
    %d1 = randperm(size(dico,1), 1);
    %d2 = randperm(size(dico,1), 1);
    %d1 = d1/norm(d1);
    %d2 = d2/norm(d2);
    %P = PMat{1}{2};
    %[u1, u2, c] = fix_angle(d1,d2,P);
    %printf('\nDoen')
    %c
    %pause
    %pause

    subdictsize  = size(dico, 2)/ nclass;
    
    % compute atom weights
    %     if nargin > 2
    %         rank = sum(amp.*amp, 2);
    %     else
    %         rank = randperm(length(dico));
    %     end
    
    % several decorrelation iterations might be needed to reach global
    % coherence mu. niter can be adjusted to needs.
    niter = 1;
    %mu = 0.5;
    %iter = 1;
    while max(max(abs((dico'*dico -eye(size(dico,2))).* B))) > mu + 0.01
        printf('\nOptimize iter %d', niter);
        % find pairs of high correlation atoms
        colors = dico_color(dico, mu, B);
        %max(max(colors))
        
        % iterate on all pairs
        nbColors = max(colors);
        isupdate = false;
        for c = 1:nbColors
            index = find(colors==c);
            if numel(index) == 2
                i1 = index(1);
                i2 = index(2);
                idx1 = ceil(index(1)/subdictsize);
                idx2 = ceil(index(2)/subdictsize);
                P = PMat{idx1}{idx2};
                u1 = dico(:,i1);
                u2 = dico(:,i2);
                %i1, i2
                c1 = centroids{idx1};
                c2 = centroids{idx2};
                [d1,d2] = fix_angle(u1, u2, P, mu, c1, c2);
                dico(:,i1) = d1;
                dico(:,i2) = d2;
                %d1'*d2
                %norm(d1) == NaN
                %printf('\n norm1 = %f, norm2 = %f, cor =%f', norm(d1), norm(d2), abs(d1'*d2))
                %if isnan(d1)
                    %printf('\nErrr')
                %end
            end
            %pause
        end
        niter = niter+1;

        %if iter 
        %imshow((dico'*dico -eye(size(dico,2))).* B);
        %pause()
    end
end

