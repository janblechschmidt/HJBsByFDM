function plotDecisionBounday(decBnd, varargin)

n = size(decBnd,1)/2;
for i = 1:n
	idx = [2*i-1,2*i];
	plot(decBnd(idx,1), decBnd(idx,2), varargin{:});
end % for i = 1:n
