function [idx, vals] = stencil_SecondOrder(a, j, d)
% Standard second order stencil for derivatives d_xx, d_yy, d_zz etc.
% INPUT: a - algorithmic structure, necessary components are spatial differences and grid sizes
%        j - vector of global indices where to apply difference operator
%        d - (d)imension controls the dimension in which the stencil works,
%            i.e. d = 1 means d_xx, d = 2 means d_yy, etc.
% OUTPUT: idx  - offset for assembling of differential operator
%         vals - coefficients of size = (ndofs - 1) x 3

if nargin < 3 || isempty(d)
	d = 1;
end % if nargin < 3 || isempty(d)

if nargin < 2 || isempty(j)
	switch d
	case 1
		j = a.ixinner;
	case 2
		j = a.iyinner;
	case 3
		j = a.izinner;
	end % switch d
end % if nargin < 2 || isempty(j)

if size(a.XY,2)==1
    h = a.Dx;
    fak = 1;
else
    switch d
    case 1
	    h = a.Dx;
    	fak = (a.ny+1);
    case 2
    	h = a.Dy;
	    fak = 1;
    case 3
	    h = a.Dz;
    	fak = (a.nx+1) * (a.ny+1);
    end % switch d
end % if size(a.XY,2)==1

% Prepare indices
jb = j - fak;
jf = j + fak;
idx = [jb, j, jf];

% Prepare coefficient values
valm = +2 ./ ( h(jb) .* ( h(jb) + h(j) ) );
val0 = -2 ./ ( h(jb) .* h(j) );
valp = +2 ./ ( h(j)  .* ( h(jb) + h(j) ) );
vals = [valm, val0, valp];

end % function
