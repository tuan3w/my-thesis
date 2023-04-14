function s = unbuffer(obj)

%% Check inputs and Defaults
if ~isprop(obj,'bufferOperator') || isempty(obj.bufferOperator)
	error('You must buffer a signal before unbuffer it, come on!');
end

switch lower(obj.bufferOperator.method)
	%Unbuffer using overlap-add method
	case 'standard'
		w = obj.bufferOperator.window(obj.bufferOperator.wLength);
		S = diag(1./w)*obj.S;
		%Non overlapping case
		if obj.bufferOperator.overlap == 0
			s = S(:);
		%Overlapping case
		else
			Stemp = S(1:obj.bufferOperator.wLength-obj.bufferOperator.overlap,1:end);
			s = [Stemp(:); S(obj.bufferOperator.wLength-obj.bufferOperator.overlap+1:end,end)];
		end
	%Unbuffer using inverse lot with identity local transform
	case 'lot'
		s = ilot(obj.S(:),obj.bufferOperator.wLength,'id',...
			obj.bufferOperator.overlap,obj.bufferOperator.window);
	otherwise
		error('Please specify a valid buffer method');
end
