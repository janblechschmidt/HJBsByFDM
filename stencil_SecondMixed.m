function [idx, vals] = stencil_SecondMixed(a, j, mode, d1, d2)
% Second order stencil for mixed derivatives of order 2
% INPUT: a - algorithmic structure, necessary components are spatial differences and grid sizes
%        j - vector of global indices where to apply difference operator
%        mode - 'Forward'  (coefficient a_12 > 0) or
%               'Backward' (coefficient a_12 < 0)
%        d1 - controls the 1st dimension in which the stencil works
%        d2 - controls the 2nd dimension in which the stencil works
% OUTPUT: idx  - offset for assembling of differential operator
%         vals - coefficients of size = (ndofs - 1) x 3

% TODO: Incorporate ordinary second order stencil

if d1 == d2
	error('Use the function stencil_SecondOrder instead.')
end % if d1 == d2

switch d1
case 1
	h1 = a.Dx;
	fak1 = (a.ny+1);
case 2
	h1 = a.Dy;
	fak1 = 1;
case 3
	h1 = a.Dz;
	fak1 = (a.nx+1) * (a.ny+1);
end % switch d1

switch d2
case 1
	h2 = a.Dx;
	fak2 = (a.ny+1);
case 2
	h2 = a.Dy;
	fak2 = 1;
case 3
	h2 = a.Dz;
	fak2 = (a.nx+1) * (a.ny+1);
end % switch d2

% Prepare indices
jleft  = j - fak1;
jright = j + fak1;
jdown  = j - fak2;
jup    = j + fak2;
jupleft    = jup - fak1;
jupright   = jup + fak1;
jdownleft  = jdown - fak1;
jdownright = jdown + fak1;

% Prepare factors
alpha = 1./( h1(j) .* h2(j) );
beta  = 1./( h1(jleft) .* h2(jdown) );
gamma = 1./( h1(j) .* h2(jdown) );
delta = 1./( h1(jleft) .* h2(j) );
D_1pD_2p.idx = [j, jright, jup, jupright];
D_1mD_2m.idx = [j, jleft, jdown, jdownleft];
D_1pD_2m.idx = [j, jright, jdown, jdownright];
D_1mD_2p.idx = [j, jleft, jup, jupleft];

D_1pD_2p.vals = [alpha, -alpha, -alpha, alpha];
D_1mD_2m.vals = [beta, -beta, -beta, beta];
D_1pD_2m.vals = [-gamma, gamma, gamma, -gamma];
D_1mD_2p.vals = [-delta, delta, delta, -delta];
% TODO Check this, there was an error in stencil_SecondOrder.m
switch mode

	case 'Forward' % (coefficient a_12 > 0)
		idx = [D_1pD_2p.idx, D_1mD_2m.idx];
		vals = 0.5*[D_1pD_2p.vals, D_1mD_2m.vals];

	case 'Backward' % (coefficient a_12 < 0)
		idx = [D_1pD_2m.idx, D_1mD_2p.idx];
		vals = 0.5*[D_1pD_2m.vals, D_1mD_2p.vals];
		% idx = [D_1pD_2m.idx, D_1pD_2m.idx];
		% vals = 0.5*[D_1pD_2m.vals, D_1pD_2m.vals];

end % function

