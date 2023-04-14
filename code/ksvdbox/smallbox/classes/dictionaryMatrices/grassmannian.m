function [A G res muMin] = grassmannian(n,m,nIter,dd1,dd2,initA,verb)
% grassmanian attempts to create an n by m matrix with minimal mutual
% coherence using an iterative projection method.
%
% [A G res] = grassmanian(n,m,nIter,dd1,dd2,initA)
%
% REFERENCE
% M. Elad, Sparse and Redundant Representations, Springer 2010.

%% Parameters and Defaults
error(nargchk(2,7,nargin));

if ~exist('verb','var')  || isempty(verb),  verb = false; end %verbose output
if ~exist('initA','var') || isempty(initA), initA = randn(n,m); end %initial matrix
if ~exist('dd2','var')   || isempty(dd2),   dd2 = 0.9; end %shrinking factor
if ~exist('dd1','var')   || isempty(dd1),   dd1 = 0.9; end %percentage of coherences to be shrinked
if ~exist('nIter','var') || isempty(nIter), nIter = 10; end %number of iterations

%% Main algo
A = normc(initA); %normalise columns
[Uinit Sigma] = svd(A);
G = A'*A; %gram matrix

muMin = sqrt((m-n)/(n*(m-1)));              %Lower bound on mutual coherence (equiangular tight frame)
res = zeros(nIter,1);
if verb
	fprintf(1,'Iter		mu_min		mu \n');
end

% optimise gram matrix
for iIter = 1:nIter
	gg  = sort(abs(G(:))); %sort inner products from less to most correlated
	pos = find(abs(G(:))>=gg(round(dd1*(m^2-m))) & abs(G(:)-1)>1e-6); %find large elements of gram matrix
	G(pos) = G(pos)*dd2;	%shrink large elements of gram matrix
	[U S V] = svd(G);	%compute new SVD of gram matrix
	S(n+1:end,1+n:end) = 0; %set small eigenvalues to zero (this ensures rank(G)<=d)
	G = U*S*V';			%update gram matrix
	G = diag(1./abs(sqrt(diag(G))))*G*diag(1./abs(sqrt(diag(G)))); %normalise gram matrix diagonal
	if verb
		Geye = G - eye(size(G));
		fprintf(1,'%6i    %12.8f  %12.8f  \n',iIter,muMin,max(abs(Geye(:))));
	end
end

[V_gram Sigma_gram] = svd(G);				%calculate svd decomposition of gramian
Sigma_new = sqrt(Sigma_gram(1:n,:)).*sign(Sigma); %calculate singular values of dictionary
A = Uinit*Sigma_new*V_gram';				%update dictionary
