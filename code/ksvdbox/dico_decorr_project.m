function dico = dico_decorr_symetric(dico, mu, B)
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
    theta = acos(mu)/2;
    ctheta = cos(theta);
    stheta = sin(theta);
    
    % compute atom weights
    %     if nargin > 2
    %         rank = sum(amp.*amp, 2);
    %     else
    %         rank = randperm(length(dico));
    %     end
    
    % several decorrelation iterations might be needed to reach global
    % coherence mu. niter can be adjusted to needs.
    niter = 1;
    %size(dico)
    %size(dictlabel)
    %size(dico)
    %printf('\nfsfjjjjjjjjjjjjj....................')
    %imshow(dico'*dico -eye(size(dico,2)));
    %pause()
    %imshow((dico'*dico -eye(size(dico,2))).* B);
    %pause()

    %iter = 1;
    while max(max(abs((dico'*dico -eye(size(dico,2))).* B))) > mu + 0.01
        %printf('\nOptimize iter %d', niter);
        % find pairs of high correlation atoms
        colors = dico_color(dico, mu, B);
        
        % iterate on all pairs
        nbColors = max(colors);
        for c = 1:nbColors
            index = find(colors==c);
            if numel(index) == 2
                if dico(:,index(1))'*dico(:,index(2)) > 0               
                    %build the basis vectors
                    v1 = dico(:,index(1))+dico(:,index(2));
                    v1 = v1/norm(v1);
                    v2 = dico(:,index(1))-dico(:,index(2));
                    v2 = v2/norm(v2);
                    
                    dico(:,index(1)) = ctheta*v1+stheta*v2;
                    dico(:,index(2)) = ctheta*v1-stheta*v2;
                else
                    v1 = dico(:,index(1))-dico(:,index(2));
                    v1 = v1/norm(v1);
                    v2 = dico(:,index(1))+dico(:,index(2));
                    v2 = v2/norm(v2);
                    
                    dico(:,index(1)) = ctheta*v1+stheta*v2;
                    dico(:,index(2)) = -ctheta*v1+stheta*v2;
                end
            end
        end
        niter = niter+1;

        %if iter 
        %imshow((dico'*dico -eye(size(dico,2))).* B);
        %pause()
    end
end

