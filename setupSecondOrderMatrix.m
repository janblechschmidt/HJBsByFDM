function A = setupSecondOrderMatrix(p, a)
% setupSecondOrderMatrix sets up the matrix for all second-order derivative terms
% defined in a problem structure p and an algorithmic structure a

% Load helper functions
helperFunctions;

% Initialize second order matrix
m = p.ndofs;
n = p.ndofs;
A = sparse(m, n);

% Evaluate diffusion coeffcients at grid points
vals = cellfun(@(f) f(a.XY), p.diffusion,'UniformOutput',false);
switch p.dim
case 1

	% D_xx derivative
	if any(vals{1,1} ~= 0)
		A2xx = setupDifferentialOperator(p, a, 'SecondOrder', 1);
		A = A + spdiags(vals{1,1},0,m,n) * A2xx;
	end

case 2

	if ~all(vals{1,2}-vals{2,1} == 0)
		error('Diffusion matrix is not symmetric\n.')
    end

	if ~all(vals{1,1}>abs(vals{1,2}))
		warning('Diffusion matrix is not diagonally dominant\n.')
	end

	% D_xx derivative
	if any(vals{1,1} ~= 0)
		A2xx = setupDifferentialOperator(p, a, 'SecondOrder', 1);
		A = A + spdiags(vals{1,1},0,m,n) * A2xx;
	end

	% D_yy derivative
	if any(vals{2,2} ~= 0)
		A2yy = setupDifferentialOperator(p, a, 'SecondOrder', 2);
		A = A + spdiags(vals{2,2},0,m,n) * A2yy;
    end

	% D_xy derivative
	if any(vals{1,2} ~=0)
		A2xyForward  = setupDifferentialOperator(p, a, 'Mixed', 1, 2, 'Forward');
		A2xyBackward = setupDifferentialOperator(p, a, 'Mixed', 1, 2, 'Backward');
        A = A + 2 * spdiags(posPart(vals{1,2}), 0, m, n) * A2xyForward;
        A = A - 2 * spdiags(negPart(vals{1,2}), 0, m, n) * A2xyBackward;
        % A = A + 2 * spdiags(posPart(vals{1,2}), 0, m, n) * A2xyForward;
        % A = A - 2 * spdiags(negPart(vals{1,2}), 0, m, n) * A2xyBackward;
	end

case 3

	if ~all(vals{1,2}-vals{2,1} == 0)
		error('Diffusion matrix is not symmetric in x and y\n.')
	end

	if ~all(vals{1,3}-vals{3,1} == 0)
		error('Diffusion matrix is not symmetric in x and z\n.')
	end

	if ~all(vals{2,3}-vals{3,2} == 0)
		error('Diffusion matrix is not symmetric in y and z\n.')
    end
    c1 = ~all(vals{1,1}>abs(vals{1,2})+abs(vals{1,3}));
    c2 = ~all(vals{2,2}>abs(vals{2,1})+abs(vals{2,3}));
    c3 = ~all(vals{3,3}>abs(vals{3,1})+abs(vals{3,2}));
	if c1 || c2 || c3
		warning('Diffusion matrix is not diagonally dominant\n.')
	end

	% D_xx derivative
	if any(vals{1,1} ~= 0)
		A2xx = setupDifferentialOperator(p, a, 'SecondOrder', 1);
		A = A + spdiags(vals{1,1},0,m,n) * A2xx;
	end

	% D_yy derivative
	if any(vals{2,2} ~= 0)
		A2yy = setupDifferentialOperator(p, a, 'SecondOrder', 2);
		A = A + spdiags(vals{2,2},0,m,n) * A2yy;
	end

	% D_zz derivative
	if any(vals{3,3} ~= 0)
		A2zz = setupDifferentialOperator(p, a, 'SecondOrder', 3);
		A = A + spdiags(vals{3,3},0,m,n) * A2zz;
	end

	% D_xy derivative
	if any(vals{1,2} ~=0)
		A2xyForward  = setupDifferentialOperator(p, a, 'Mixed', 1, 2, 'Forward');
		A2xyBackward = setupDifferentialOperator(p, a, 'Mixed', 1, 2, 'Backward');
		A = A + 2 * spdiags(posPart(vals{1,2}), 0, m, n) * A2xyForward;
		A = A - 2 * spdiags(negPart(vals{1,2}), 0, m, n) * A2xyBackward;
	end

	% D_xz derivative
	if any(vals{1,3} ~=0)
		A2xzForward  = setupDifferentialOperator(p, a, 'Mixed', 1, 3, 'Forward');
		A2xzBackward = setupDifferentialOperator(p, a, 'Mixed', 1, 3, 'Backward');
		A = A + 2 * spdiags(posPart(vals{1,3}), 0, m, n) * A2xzForward;
		A = A - 2 * spdiags(negPart(vals{1,3}), 0, m, n) * A2xzBackward;
	end

	% D_yz derivative
	if any(vals{2,3} ~=0)
		A2xzForward  = setupDifferentialOperator(p, a, 'Mixed', 2, 3, 'Forward');
		A2xzBackward = setupDifferentialOperator(p, a, 'Mixed', 2, 3, 'Backward');
		A = A + 2 * spdiags(posPart(vals{2,3}), 0, m, n) * A2xzForward;
		A = A - 2 * spdiags(negPart(vals{2,3}), 0, m, n) * A2xzBackward;
	end

end % switch p.dim

end % function setupSecondOrderMatrix(p, a)
