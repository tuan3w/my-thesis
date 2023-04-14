function [dic mus srs] = iterativeprojections(dic,mu,Y,X,params)
% grassmanian attempts to create an n by m matrix with minimal mutual
% coherence using an iterative projection method.
%
% REFERENCE
%

%% Parameters and Defaults
if ~nargin, testiterativeprojections; return; end

if ~exist('params','var') || isempty(param), params = struct; end
if ~isfield(params,'nIter'), params.nIter = 10; end %number of iterations
if ~isfield(params,'eps'),	 params.eps = 1e-9;	 end %tolerance level
[n m] = size(dic);

SNR = @(dic) snr(Y,dic*X);									 %SNR function
MU  = @(dic) max(max(abs((dic'*dic)-diag(diag(dic'*dic))))); %coherence function

%% Main algorithm
dic = normc(dic);	%normalise columns
alpha = m/n;		%ratio between number of atoms and ambient dimension

mus = zeros(params.nIter,1);		%coherence at each iteration
srs = zeros(params.nIter,1);		%signal to noise ratio at each iteration
iIter = 1;
while iIter<=params.nIter && MU(dic)>mu
	fprintf(1,'Iteration number %u\n', iIter);
	% calculate snr and coherence
	mus(iIter) = MU(dic);
	srs(iIter) = SNR(dic);
	
	% calculate gram matrix
	G = dic'*dic;
	
	% project onto the structural constraint set
	H = zeros(size(G));				%initialise matrix
	ind1 = find(abs(G(:))<=mu);		%find elements smaller than mu
	ind2 = find(abs(G(:))>mu);		%find elements bigger than mu
	H(ind1) = G(ind1);				%copy elements belonging to ind1
	H(ind2) = mu*sign(G(ind2));		%threshold elements belonging to ind2
	H(1:m+1:end) = 1;				%set diagonal to one
	
	% project into spectral constraint set
	[~ , S, V] = svd(H);
	%G = alpha*(V(:,1:n)*V(:,1:n)');
	G = V(:,1:n)*S(1:n,1:n)*V(:,1:n)';
	
	% calculate dictionary
	[~, S V] = svd(G);
	dic = sqrt(S(1:n,:))*V';
	
	% rotate dictionary
	options = struct('nIter',100,'step',0.001);
	[~, ~, W] = rotatematrix(Y,dic*X,'conjgradlie',options);
	dic = W*dic;
	
	iIter = iIter+1;
end
if iIter<params.nIter
	mus(iIter:end) = mus(iIter);
	srs(iIter:end) = srs(iIter);
end
% Test function
function testiterativeprojections
clc
%define parameters
n = 256;							%ambient dimension
m = 512;							%number of atoms
N = 1024;							%number of signals
mu_min = sqrt((m-n)/(n*(m-1)));		%minimum coherence

%initialise data
X = sprandn(m,N,1);					%matrix of coefficients
phi = normc(randn(n,m));			%dictionary
temp = randn(n);
W = expm(0.5*(temp-temp'));			%rotation matrix
Y = W*phi*X;						%observed signals

%optimise dictionary
[~, mus srs] = iterativeprojections(phi,0.2,Y,X);

%plot results
nIter = length(mus);

figure, 
subplot(2,1,1)
plot(1:nIter,srs,'kd-');
xlabel('nIter');
ylabel('snr (dB)');
grid on

subplot(2,1,2), hold on
plot(1:nIter,mus,'ko-');
plot([1 nIter],[mu_min mu_min],'k')
grid on
legend('\mu','\mu_{min}');

