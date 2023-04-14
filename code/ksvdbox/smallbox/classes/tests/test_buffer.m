%% Buffer Test
% Test script for the function audio->buffer which takes an audio object as
% an input and buffers it into frames. The test assesses wether all the
% possible buffer methods are invertible.
%%

N = 500;					% number of audio samples
TOL = 200;					% tolerance (in decibels)
verb = true;				% verbose

obj = audio(randn(N,1));	% audio object
% wLengts: one window, two windows, maximum number of windows
temp = factor(N);
wLengths = [N; fix(N/2); fix(N/temp(end))];
% overlaps: zero, half window, max
overlapNames = {'zero','half-window','maximum'};
overlaps = {@(n)0, @(n)n/2, @(n)n-1};
% windows: valid windows
windows = {@hanning,@hamming,@triang,@rectwin};
% methods: valid methods
methods = {'standard'};


nLen  = length(wLengths);
nOver = length(overlaps);
nWin  = length(windows);
nMet  = length(methods);
count = 1;
for iLen=1:nLen
	for iOver=1:nOver
		for iWin=1:nWin
			for iMet=1:nMet
				if verb
				printf('\n buffer test %d/%d - %d window length, %s overlap, %s window and %s method ... ',...
					count,nMet*nWin*nOver*nLen,wLengths(iLen),overlapNames{iOver},func2str(windows{iWin}),methods{iMet});
				end
				obj = buffer(obj,wLengths(iLen),overlaps{iOver}(wLengths(iLen)),windows{iWin},methods{iMet});
				s_rec = obj.unbuffer;
				if snr(obj.s,s_rec) > TOL
					if verb, printf('Passed'); end
				else
					error('Test failed!');
				end
				count = count+1;
			end
		end
	end
end
