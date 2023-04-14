function mu = cumcoherence(obj,p)
% Calculates the p-cumulative coherence of the dictionary, defined as
% \mu_p(k) = max_|I|=k{max_j\notin I{(\sum_i \in I}|<\phi_i,\phi_j>|^p)^1/p}
%
% INPUT
% obj: dictionary object
% p: power (default 1)
%
% OUTPUT
% mu: p-cumulative coherence
if ~exist('p','var') || isempty(p), p = 1; end

obj = normalize(obj);
[M N] = size(obj.phi);
mu = zeros(M,1);
for m=1:M
	c = zeros(N);
	for i=1:N
		c(:,i) = abs(obj.phi'*obj.phi(:,i)).^p;
		c(i,i) = 0;
	end
	c = sort(c,'descend');
	c = c(1:m,:);
	if m==1
		mu(m) = max(c.^(1/p));
	else
		mu(m) = max(sum(c).^(1/p));
	end
end
end
