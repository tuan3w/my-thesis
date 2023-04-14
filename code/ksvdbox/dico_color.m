function colors = dico_color(dico, mu, Mask)
    % DICO_COLOR cluster a dictionary in pairs of high correlation atoms.
    % Called by dico_decorr.
    %
    % Parameters:
    % -dico: the dictionary
    % -mu: the correlation threshold
    %
    % Result:
    % -colors: a vector of indices. Two atoms with the same color have a 
    % correlation greater than mu 
    
    numAtoms = length(dico);
    colors = zeros(numAtoms, 1);
    
    % compute the correlations
    G = abs(dico'*dico);
    G = G-eye(size(G));
    G = G.* Mask;
    
    % iterate on the correlations higher than mu
    c = 1;   
    maxCorr = max(max(G));
    while maxCorr > mu
        % find the highest correlated pair
        x = find(max(G)==maxCorr, 1);
        y = find(G(x,:)==maxCorr, 1);
        
        % color them if they are from difference class
        %if label(x) ~= label(y)
        colors(x) = c;
        colors(y) = c;
        c = c+1;
        %end
        
        % make sure these atoms never get selected again
        G(x,:) = 0;
        G(:,x) = 0;
        G(y,:) = 0;
        G(:,y) = 0;
        
        % find the next correlation
        maxCorr = max(max(G));
    end
    
    % complete the coloring with singletons
    index = find(colors==0);
    colors(index) = c:c+length(index)-1;
end
