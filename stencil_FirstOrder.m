function [idx, vals] = stencil_FirstOrder(a, mode, j, d)
% Collection of stencils for first order approximation
% INPUT: a    - (a)lgorithmic structure, necessary components are spatial differences and grid sizes
%        mode - option to choose between 'Forward', 'Backward', 'Central'
%        j    - vector of global indices where to apply difference operator
%        d    - (d)imension controls the dimension in which the stencil works,
%               i.e. d = 1 means d_xx, d = 2 means d_yy, etc.
% OUTPUT: idx  - offset for assembling of differential operator
%         vals - coefficients of size = (ndofs - 1) x 3

if nargin < 4 || isempty(d)
	d = 1;
end % if nargin < 4

if nargin < 3 || isempty(j)
	switch mode
	case 'Forward'
		switch d
		case 1
			j = setdiff(a.iall, a.ixup);
		case 2
			j = setdiff(a.iall, a.iyup);
		case 3
			j = setdiff(a.iall, a.izup);
		end % switch d

	case 'Backward'
		switch d
		case 1
			j = setdiff(a.iall, a.ixlow);
		case 2
			j = setdiff(a.iall, a.iylow);
		case 3
			j = setdiff(a.iall, a.izlow);
		end % switch d

	case 'Central'
		switch d
		case 1
			j = setdiff(a.iall, unique([a.ixlow; a.ixup]));
		case 2
			j = setdiff(a.iall, unique([a.iylow; a.iyup]));
		case 3
			j = setdiff(a.iall, unique([a.izlow; a.izup]));
		end % switch d
	end % switch mode
end % if nargin < 3 || isempty(j)


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

% Forward first order stencil which is of order 1
if strcmp(mode,'Forward')
	idx = [j, jf];
	vals = 1./h(j) * [-1, 1];
end

% Backward first order stencil which is of order 1
if strcmp(mode,'Backward')
	idx = [jb, j];
	vals = 1./h(jb) * [-1, 1];
end

% Central first order stencil which is of order 2
if strcmp(mode,'Central')
	idx = [jb, jf];
	vals = 1./( h(jb) + h(j) ) * [-1, 1];
end
end % function
