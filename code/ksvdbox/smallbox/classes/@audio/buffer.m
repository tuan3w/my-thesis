%% Buffer function
% Buffers the samples of the audio object into the columns of the matrix S
% based on the input parameters
%%
function obj = buffer(obj,wLength,overlap,window,method)

%% Check inputs & Defaults
error(nargchk(2, 5, nargin, 'struct'));
if rem(length(obj.s),wLength)
% 	error('The wLength must be an integer divisor of the signal length!'); 
end
if ~exist('overlap','var') || isempty(overlap), overlap = 0; end
if ~exist('method','var') || isempty(method), method = 'standard'; end

%% Buffer audio
switch lower(method)
	case 'standard'
		if ~exist('window','var') || isempty(window), window = @rectwin; end
		validWindows = {'hanning','hamming','triang','rectwin'};
		if ~sum(strcmpi(validWindows,func2str(window)));
			error('The window chosen is not valid because it cannot be inverted!');
		end
		obj.S = diag(window(wLength))*buffer(obj.s,wLength,overlap,'nodelay');
% 	case 'lot'
% 		if ~exist('window','var') || isempty(window), window = 'sin2'; end
% 		s_lot = lot(obj.s,wLength,'id',overlap,window);
% 		obj.S = buffer(s_lot,wLength);
	otherwise
		error('Please specify a valid buffer method');
end

obj.bufferOperator = struct('wLength',wLength,'overlap',...
	overlap,'window',window,'method',method);
