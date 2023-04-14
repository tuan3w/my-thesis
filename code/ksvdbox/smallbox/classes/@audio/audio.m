%% AUDIO OBJECT CLASS
% Class designed to analyse and process audio signals
%%
classdef audio
	properties (SetAccess = protected)
		s				%vector containing the audio signal
		fs				%sampling frequency
		nBits			%number of bits per sample
		name			%string containing the name of the audio file
		format			%string containing the format of the audio file
		bufferOperator	%struct containing the parameters of the buffer operator
		S				%matrix containing frames of audio
	end
	
	methods
		%% Constructor
		function obj = audio(varargin)
			%%% obj = audio(varargin)
			% Audio object constructor.
			% INPUT: either a path to an audio file, or the following
			% arguments.
			% - s: vector containing the audio samples
			% - fs: sampling frequency
			% - nBits: number of bits per sample
			% - name: name of the audio object
			% - format: format of the audio object
			%
			% if no arguments are specified, prompt for the choice of an
			% audio file
			if ~nargin
				[fileName,pathname] = uigetfile({'*.wav; *.aiff;'},'Select an audio file');
				varargin{1} = strcat(pathname,filesep,fileName);
			end
			% if a file is specified, read it from disk
			if ischar(varargin{1})
				[~, obj.name obj.format] = fileparts(varargin{1});
				switch obj.format
					case '.wav'
						[obj.s obj.fs obj.nBits] = wavread(varargin{1});
					otherwise
						error('Unsupported audio format')
				end
			% if properties are specified, set them to input values
			else
				obj.s = varargin{1};
				if nargin>1, obj.fs = varargin{2}; else obj.fs = []; end
				if nargin>2, obj.nBits = varargin{3}; else obj.nBits = []; end
				if nargin>3, obj.name = varargin{4}; else obj.name = []; end
				if nargin>4, obj.format = varargin{5}; else obj.format = []; end
			end
			obj.S = [];
			obj.bufferOperator = [];
		end
		
		%% Playback functions
		function player = play(obj, player)
			if ~exist('player','var') || isempty(player)
				player = audioplayer(obj.s,obj.fs);
			end
			play(player);
		end
		
		function player = stop(obj, player)
			if ~exist('player','var') || isempty(player)
				player = audioplayer(obj.s,obj.fs);
			end
			stop(player)
		end
		
		function player = pause(obj, player)
			if ~exist('player','var') || isempty(player)
				player = audioplayer(obj.s,obj.fs);
			end
			pause(player)
		end
	end
end
