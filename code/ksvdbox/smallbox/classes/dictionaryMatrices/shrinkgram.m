function [dic mus] = shrinkgram(dic,mu,dd1,dd2,params)
% grassmanian attempts to create an n by m matrix with minimal mutual
% coherence using an iterative projection method.
%
% [A G res] = grassmanian(n,m,nIter,dd1,dd2,initA)
%
% REFERENCE
% M. Elad, Sparse and Redundant Representations, Springer 2010.

%% Parameters and Defaults
if ~nargin, testshrinkgram; return; end

if ~exist('dd2','var')   || isempty(dd2),   dd2 = 0.9; end %shrinking factor
if ~exist('dd1','var')   || isempty(dd1),   dd1 = 0.9; end %percentage of coherences to be shrinked
if ~exist('params','var') || isempty(params), params = struct; end
if ~isfield(params,'nIter'), params.nIter = 100; end

%% Main algo
dic = normc(dic);	%normalise columns
G = dic'*dic;		%gram matrix
[n m] = size(dic);

MU  = @(G) max(max(abs(G-diag(diag(G))))); %coherence function

mus   = ones(params.nIter,1);
iIter = 1;
% optimise gram matrix
while iIter<=params.nIter && MU(G)>mu
	mus(iIter) = MU(G);		%calculate coherence
	gg  = sort(abs(G(:)));	%sort inner products from less to most correlated
	pos = find(abs(G(:))>=gg(round(dd1*(m^2-m))) & abs(G(:)-1)>1e-6); %find large elements of gram matrix
	G(pos) = G(pos)*dd2;	%shrink large elements of gram matrix
	[U S V] = svd(G);		%compute new SVD of gram matrix
	S(n+1:end,1+n:end) = 0; %set small eigenvalues to zero (this ensures rank(G)<=d)
	G = U*S*V';				%update gram matrix
	G = diag(1./abs(sqrt(diag(G))))*G*diag(1./abs(sqrt(diag(G)))); %normalise gram matrix diagonal
	iIter = iIter+1;
end
%if iIter<params.nIter
%	mus(iIter:end) = mus(iIter-1);
%end

[V_gram Sigma_gram] = svd(G);				%calculate svd decomposition of gramian
dic = sqrt(Sigma_gram(1:n,:))*V_gram';		%update dictionary

function testshrinkgram
clc
%define parameters
n = 256;							%ambient dimension
m = 512;							%number of atoms
N = 1024;							%number of signals
mu_min = sqrt((m-n)/(n*(m-1)));		%minimum coherence

%initialise data
phi = normc(randn(n,m));			%dictionary

%optimise dictionary
[~, mus] = shrinkgram(phi,0.2);

%plot results
nIter = length(mus);

figure, hold on
plot(1:nIter,mus,'ko-');
plot([1 nIter],[mu_min mu_min],'k')
grid on
legend('\mu','\mu_{min}');
