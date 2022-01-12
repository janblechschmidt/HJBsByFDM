function [x,nx] = refineCells(x, indices, cnt)
% refineCells refines the cells from x(indices(1)) to x(indices(end))
if nargin < 3
	cnt = 2;
end

xnew = x;
for i=1:(cnt-1)
	xtmp = x(indices(1:end-1)) + i/cnt * (x(indices(2:end)) - x(indices(1:end-1)) );
	xnew = union(xnew,xtmp);
end

x = xnew;
nx = length(x) - 1;
