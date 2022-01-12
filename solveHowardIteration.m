function [v, u] = solveHowardIteration(p, a, v_prev, u_prev)

% Prepare Howard loop
done = 0;
iter = 0;

check_inverse_monotonicity = 0;

w_prev = v_prev;
w_linferror_prev = inf;
w_l2error_prev = inf;

Id = speye(p.ndofs);

% Print head of convergence table on each 10th iteration
if mod(p.global_iter, 10) == 0
	fprintf('\n   Iter    |   Time   | ||w_n-w_nm1||_inf   |   ||w_n - w_nm1||_2  |   ||u_n - u_nm1||_2\n')
    if p.trivialControl
        fprintf('----------------------------------------------------------------------------------------------\n')
    end
end
if ~p.trivialControl
    fprintf('----------------------------------------------------------------------------------------------\n')
end


switchingPoints = [];
while done == 0

    % Determine optimal control and matrix for first order term using this control
    if p.trivialControl
        u = zeros(p.ndofs,p.control_dim);
        A1 = a.A1func(p.t, u);
        done = 1;
        xi = 0;
    else
	    [u, ~, A1, xi] = p.determineOptimalControl(p, a, w_prev, v_prev);
    end
    
	% Set up source term for current control u
    F = p.f(a.XY, u);

	% Setup system matrix
	S = Id - p.theta * p.ht * (a.A2 + A1 - a.A0);

    if check_inverse_monotonicity
        S_dof = S(a.dof_idxs, a.dof_idxs);
        % s = 0;
        % for i =1:length(a.dof_idxs)
        %     ei = zeros(length(a.dof_idxs),1);
        %     ei(i) = 1;
        %     s = s+ min(S_dof\ei)<-1e-15;
        % end
        % s
        Sinv = inv(S_dof);
    
        smin = min(min(Sinv));
        smax = max(max(Sinv));
        if  smin < -1e-15
            fprintf('S in not inverse monotone [%4.2f, %4.2f]\n', full(smin), full(smax));
            keyboard
        end
    end

	% Setup rhs
	rhs = (Id + (1 - p.theta) * p.ht * (a.A2 + A1 - a.A0) ) * v_prev + p.ht * F;

    % Incorporate active set, e.g., in American option pricing setting
    if xi == 1
        active_idxs = find(u==1);
        S(active_idxs,:) = 0;
		S = S + sparse(active_idxs, active_idxs, ones(size(active_idxs)), p.ndofs, p.ndofs);
        rhs(active_idxs) = p.boundaryVal(p.Tmax, a.XY(active_idxs,:));
    end
    save('tmp.mat', 'S')
    % keyboard


	% Incorporate Dirichlet boundary conditions, if necessary
	if p.DirichletBC

		% Get boundary values
		bndval = p.boundaryVal(p.t, a.XY(a.bc_idxs, :));

		% Incorporate Dirichlet boundary conditions
		S(a.bc_idxs,:) = 0;
		S = S + sparse(a.bc_idxs, a.bc_idxs, ones(size(a.bc_idxs)), p.ndofs, p.ndofs);
		rhs(a.bc_idxs) = bndval;

	end % if p.DirichletBC

	if ~isStrictlyDiagonallyDominant(S)
		fprintf('Matrix S is not strictly diagonally dominant\n')
        keyboard
	else
		% fprintf('Matrix S is strictly diagonally dominant\n')
    end

    w = S \ rhs;
    % keyboard

	w_linferror = a.norm1(w - w_prev);
	w_l2error = a.norm2(w - w_prev);
    u_l2error = a.norm2(u - u_prev);
    qoi_val = a.get_qoi_val(w);
    fprintf(' %4d (%2d) |   %4.2f   |       %8.2e      |       %8.2e       |',p.global_iter, iter, p.t, w_linferror, w_l2error)
    for i = 1:p.control_dim
        fprintf(' %8.2e ', u_l2error(i))
    end
    fprintf('| %8.2e \n', qoi_val);
	if p.plotLevel > 2
		sfigure(1); clf;
		plotValue(p, a, w)

		sfigure(2); clf;
		plotControl(p, a, u)
        pause
	end % if p.plotLevel > 2

	v = w;
	if (w_linferror < p.errorLevelAbsolute) || (iter > p.howard_maxiter) || p.theta == 0
		done = 1;
	else
		w_prev = w;
		u_prev = u;
		w_linferror_prev = w_linferror;
		w_l2error_prev = w_l2error;
	end % if (w_linferror < p.errorLevelAbsolute) || (iter > p.howard_maxiter)


	iter = iter + 1;
end % while done == 0
