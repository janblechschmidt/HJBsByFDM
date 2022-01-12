function [v, u] = solveSemiLagrangeIteration(p, a, v_prev, u_prev)

% 0th step: Initialize explicit and implicit system matrix
% These matrices don't change anymore with varying controls,
% but can vary due to temporal changes

Id = speye(p.ndofs);

% Set up matrix for convective term in stochastic directions
A1 = a.A1func(p.t, zeros(size(u_prev)));


% Get up matrix for convective term for zero control u = 0
% switch p.dim
% case 2
%     A1 = a.A1func(p.t, zeros(size(u_prev)));
% case 3
%     A1 = a.A1func(p.t, +a.Z);
% end % switch p.dim

% Setup system matrix
Simp = Id - p.ht * p.theta * (a.A2 + A1 - a.A0);
Sexp = Id + p.ht * (1 - p.theta) * (a.A2 + A1 - a.A0);

% Print head of convergence table on each 10th iteration
if mod(p.global_iter, 10) == 0
	fprintf('\n   Iter    |   Time   |     QOI\n')
end
fprintf('-------------------------------------\n')

if p.timeConstantConvection
    version = 1;
else
    version = 2;
end

switch p.dim
case 2
	switch version
	case 1
		M = cell2mat(cellfun(@(x,y) x*v_prev+y , a.T, a.F,'UniformOutput',false));
	case 2
		[X, Y] = meshgrid(a.x, a.y);
		myReshape = @(x) reshape(x,[p.ny+1,p.nx+1]);

		M = zeros(size(a.cntrls));
		for k = 1:size(a.cntrls,2)
			ctr = p.convection{2}(a.XY, p.t, a.cntrls(:,k));
			YafterControl = Y + p.ht * myReshape(ctr);
            % keyboard
			if min(min(YafterControl)) < p.ymin || max(max(YafterControl)) > p.ymax
				% TODO: Ensure, that we always remain inside the computational domain
				warning('YafterControl outside computational domain\n')
				% keyboard
			end
			% Version with interp2 is much faster here
			% for xi = 1:a.nx+1
			% 	idx = (a.ny+1)*(xi-1) + (1:a.ny+1)';
			% 	% VafterControl(idx) = interp1(a.Y(idx),Sexp(idx,idx)*v_prev(idx), a.Y(idx) + p.ht * ctr(idx)) + p.ht * c.f(a.XY(idx,:),ctr(idx));
			% 	M(idx,k) = interp1(a.Y(idx),Sexp(idx,idx)*v_prev(idx), a.Y(idx) + p.ht * ctr(idx)) + p.ht * c.f(a.XY(idx,:),ctr(idx));
			% end
			VafterControl = interp2(X, Y, myReshape(Sexp*v_prev), X, YafterControl) + p.ht * myReshape(p.f(a.XY, a.cntrls(:,k)));
			M(:,k) = VafterControl(:);
		end

	end % switch version

case 3
	switch version
	case 1
		M = cell2mat(cellfun(@(x,y) x*v_prev+y , a.T, a.F,'UniformOutput',false));
	case 2
		[X, Y, Z] = meshgrid(a.x, a.y, a.z);
		myReshape = @(x) p.vec2Mat(x);

		M = zeros(size(a.cntrls));

		for k = 1:size(a.cntrls,2)
			ctr = p.convection{2}(a.XY, p.t, a.cntrls(:,k));
			YafterControl = Y + p.ht * myReshape(ctr);
			% ctr = a.cntrls(:,k);
			% YafterControl = Y + p.ht * myReshape(ctr) - p.ht * Z;
			if min(min(min(YafterControl))) < p.ymin || max(max(max(YafterControl))) > p.ymax
				% TODO: Ensure, that we always remain inside the computational domain
				% TODO: Check, what interp3 returns in this case, if its nan, than this should be okay
				warning('YafterControl outside computational domain\n')
				% keyboard
			end
			VafterControl = interp3(X, Y, Z, myReshape(Sexp*v_prev), X, YafterControl, Z) + p.ht * myReshape(p.f(a.XY,a.cntrls(:,k)));
			% VafterControl = interp3(X, Y, Z, myReshape(Sexp*v_prev), X, YafterControl, Z) + p.ht * myReshape(c.f(a.XY,ctr));
			M(:,k) = VafterControl(:);
		end
		% keyboard
	end % switch version

end % switch p.dim

[rhs, idx] = max(M,[],2);

% We obtain u.
u = a.cntrls(sub2ind(size(a.cntrls),(1:size(idx,1))', idx ));

% 2nd step: Solve the system and obtain v
switch p.dim
case 2
	% Direct solution does not respect the block-structure
	% v = Simp \ rhs;

	% Solution for slices for fixed values of y are faster
	v = zeros(size(v_prev));
	for yi = 1:(a.ny+1)
		idx = yi+(a.ny+1)*(0:a.nx)';
		v(idx) = Simp(idx,idx)\rhs(idx);
	end % for yi = 1:(a.ny+1)
	% TODO: Is this error negligible?
	% v1 = Simp \ rhs;
	% max(abs(v - v1))
	% pause
case 3

	if p.parallel
		v = zeros(size(v_prev));
		clear idx
		parfor yi = 1:(a.ny+1)
			% all(idx == find(yi == a.iy))
			idx{yi} = yi+(a.ny+1)*(0:(a.nx+1)*(a.nz+1)-1)';
			sol{yi} = Simp(idx{yi},idx{yi})\rhs(idx{yi});
		end % for yi = 1:(a.ny+1)
		for yi = 1:(a.ny+1)
			v(idx{yi}) = sol{yi};
		end % for yi = 1:(a.ny+1)
			% v(idx) = Simp(idx,idx)\rhs(idx);
		% keyboard
		% v = zeros(size(v_prev));
		% keyboard
		% job = createJob(sched);
		% for yi = 1:(a.ny+1)
		% 	myFunc = @(idx) Simp(idx,idx)\rhs(idx);
		% 
		% 	idx = yi+(a.ny+1)*(0:(a.nx+1)*(a.nz+1)-1)';
		% 	v(idx) = Simp(idx,idx)\rhs(idx);
		% end % for yi = 1:(a.ny+1)
		% submit(job);
		% job.state
		% waitForState(job);

	else
		% Direct solution does not respect the block-structure
        % v = Simp \ rhs;

		% Solution for slices for fixed values of y are faster
        v = zeros(size(v_prev));
        for yi = 1:(a.ny+1)
            % all(idx == find(yi == a.iy))
            idx = yi+(a.ny+1)*(0:(a.nx+1)*(a.nz+1)-1)';
            v(idx) = Simp(idx,idx)\rhs(idx);
        end % for yi = 1:(a.ny+1)

		% Direct solution does not respect the block-structure
		% TODO: Is this error negligible?
        % v1 = Simp \ rhs;
        % max(abs(v - v1)) % ~ 1e-8
        % pause
	end % if p.parallel
end
qoi_val = a.get_qoi_val(v);

fprintf(' %4d (%2d) |  %4.2f  |  %8.2e   \n', p.global_iter, 0, p.t, qoi_val);

% 
% % Create graph from adjacency matrix A
% A = Simp;
% G = digraph(A-diag(diag(A)));
% 
% % Calculate connected components
% GG = conncomp(G);
% max(GG)
% 
% [~,ix] = sort(GG);
% 
% % Reorder matrix A
% AA = A(ix,ix);
% spy(AA)
% keyboard
