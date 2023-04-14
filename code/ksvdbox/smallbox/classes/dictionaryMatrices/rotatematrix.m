function [Dhat cost W] = rotatematrix(D,Phi,method,param)
% 
%
%
% REFERENCE
% M.D. Plumbley, Geometrical Methods for Non-Negative ICA: Manifolds, Lie
% Groups and Toral Subalgebra, Neurocomputing

%% Parse inputs and set defaults
if ~nargin, testrotatematrix; return, end

if ~exist('param','var')  || isempty(param),  param  = struct; end
if ~exist('method','var') || isempty(method), method = 'conjgradLie'; end
if ~isfield(param,'nIter'), param.nIter = 100; end	%number of iterations
if ~isfield(param,'eps'), param.eps = 1e-9; end		%tolerance level
if ~isfield(param,'step'), param.step = 0.01; end

J = @(W) 0.5*norm(D-W*Phi,'fro');		%cost function

% Initialise variables
cost = zeros(param.nIter,1);			%cost at each iteration
W	 = eye(size(Phi,1));				%rotation matrix
grad = ones(size(W));					%gradient
t	 = param.step;						%step size
Gprev = 0;								%previous gradient
Hprev = 0;								%previous Lie search direction
iIter = 1;								%iteration counter

%% Main algorithm
while iIter<=param.nIter && norm(grad,'fro')>eps
	cost(iIter) = J(W);				%calculate cost
	grad = (W*Phi-D)*Phi';			%calculate gradient
	switch method
		case 'unconstrained'		% gradient descent
			eta = param.step;
			W = W - eta*grad;		% update W by steepest descent
		case 'tangent'				% self correcting tangent
			eta = param.step;
			W = W - 0.5*eta*(grad - W*grad'*W);
			[U , ~, V] = svd(W);
			W = U*V';
		case 'steepestlie'			%steepest descent in Lie algebra
			eta = param.step;
			B = 2*skew(grad*W');	% calculate gradient in Lie algebra
			W = expm(-eta*B)*W;		% update W by steepest descent
		case 'linesearchlie'		% line search in Lie algebra
			B = 2*skew(grad*W');	% calculate gradient in Lie algebra
			H = -B;					% calculate direction as negative gradient
			t = searchline(J,H,W,t);% line search in one-parameter Lie subalgebra
			W = expm(t*H)*W;		% update W by line search
		case 'conjgradlie'			% conjugate gradient in Lie algebra
			G = 2*skew(grad*W');	% calculate gradient in Lie algebra
			H = -G + polakRibiere(G,Gprev)*Hprev; %calculate conjugate gradient direction
			t = searchline(J,H,W,t);% line search in one-parameter Lie subalgebra
			W = expm(t*H)*W;		% update W by line search
			Hprev = H;				% save search direction
			Gprev = G;				% save gradient
	end
	iIter = iIter+1;				% update iteration counter
end
Dhat = W*Phi;						%rotate matrix
cost(iIter:end) = cost(iIter-1);	%zero-pad cost vector
end

%% Support functions
function gamma = polakRibiere(G,Gprev)
%Polak-Ribiere rule for conjugate direction calculation
gamma = G(:)'*(G(:)-Gprev(:))/(norm(Gprev(:))^2);
if isnan(gamma) || isinf(gamma)
	gamma = 0;
end
end

function t = searchline(J,H,W,t)
%Line search in one-parameter Lie subalgebra
t = fminsearch(@(x) J(expm(x*H)*W),t);
end

function B = skew(A)
%Skew-symmetric matrix
B = 0.5*(A - A');
end


%% Test function
function testrotatematrix
clear, clc, close all
n = 256;							%ambient dimension
m = 512;							%number of atoms
param.nIter = 300;				%number of iterations
param.step  = 0.001;			%step size
param.mu    = 0.01;			%regularization factor (for tangent method)
methods = {'tangent','linesearchlie','conjgradlie'};

Phi = randn(n,m);				%initial dictionary
Qtrue = expm(skew(randn(n)));	%rotation matrix
D = Qtrue*Phi;					%target dictionary

cost = zeros(param.nIter,length(methods));
times = zeros(param.nIter,length(methods));
for iIter=1:length(methods)
	tic
	[~, cost(:,iIter)] = rotatematrix(D,Phi,methods{iIter},param);
	times(:,iIter) = linspace(0,toc,param.nIter);
	sprintf('Method %s completed in %f seconds \n',methods{iIter},toc)
end

figure, plot(times,cost)
set(gca,'XScale','lin','Yscale','log')
legend(methods)
grid on
xlabel('time (sec)')
ylabel('J(W)')
end
