clear all, clc
% Load helper functions such as positive, negative part
helperFunctions
p.dim = 2;
p.mode = 'MinimumArrivalTime';
p.outputdir = 'output2d';
p.trivialControl = 0;

p.prefix = 'test';

p.parallel = 0;
p.refinement.level = 0;
% p.theta = 0.5;
p.theta = 1.0;
p.howard_maxiter = 50;
p.useUpwinding = 1;
p.convectionExplicit = 0;
p.damping = 0;
p.method = 'Howard';
p.xmin = -1;
p.ymin = -1;
p.nx = 200;
p.ny = 200;

% Example 1
p.alpha = 0.0;
p.beta = 0.0;
p.sigmax = 1.0;
p.sigmay = 1.0;
% FD max(v): 0.358438198130502
% FV max(v): 0.358438198130563

% Example 2
p.alpha = 0.0;
p.beta = 0.0;
p.sigmax = 0.0;
p.sigmay = 0.0;
% FD: 0.916395480166108
% FV: 0.916395480168803

% Example 3
p.alpha = 0.5;
p.beta = 0.5;
p.sigmax = 0.5;
p.sigmay = 0.5;
% FD: 0.937344595250186
% FV: 0.937344595250332


p = example_MinimumArrivalTime_Problem(p);
p.nt = 100;
p.plotLevel  = 1;
p.writeLevel = 0;
p.writeValueFunction = 0;
p.writeControlFunction = 0;
p.storeLevel = 0;
p.postprocessLevel = 0;
p.errorLevelAbsolute = 1e-8;
[a, p] = setupStructure(p);

a.norm1 = @(v) norm_Linfty(v);
a.norm2 = @(v) norm_L2(a.Hv, v);

if p.plotLevel > 3
	sfigure(3); clf
	plotDiffusion(p, a);
	sfigure(4); clf
	plotConvection(p, a, 0, zeros(size(a.X)));
	% pause
end

% Setup second order matrix
a.A2 = setupSecondOrderMatrix(p, a);

% Setup routines to set up first order matrix/matrices
[a.A1func] = setupFirstOrderMatrix(p, a);

a.A0 = spdiags(p.potential(a.XY), 0, p.ndofs, p.ndofs);

% Solve 
% -v_t - sup_{u \in U} {A(x): (nabla^2 v)(x) + B(x)' * \grad(v)(x) - C(x) + F(x)} = 0

m = p.ndofs;
n = p.ndofs;
Id = speye(m);

v_prev = p.finalTimeVal(a.XY);
u_prev = zeros(size(v_prev));

% a.Ffunc = @(u) p.f(a.XY,u);


if p.storeLevel
	U = zeros(p.ndofs, p.nt+1);
	V = zeros(p.ndofs, p.nt+1);
	% Ind = zeros(p.ndofs, p.nt+1);
	% Vy = zeros(p.ndofs, p.nt+1);
end

p.global_iter = 0;

if p.writeLevel > 1
	exportVTKSolution(p, a, v_prev, u_prev, sprintf('solution_%06d.vtk',p.global_iter) );
end % if p.writeLevel > 1

% Start temporal loop
for i = 1:p.nt

	% Set current time
	p.t = a.t(p.nt+1-i);


	if strcmp(p.method, 'SemiLagrangian')
		[v, u] = solveSemiLagrangeIteration(p, a, v_prev, u_prev);
	else
		[v, u] = solveHowardIteration(p, a, v_prev, u_prev);
	end

	if p.plotLevel > 1
		sfigure(1); clf;
		plotValue(p, a, v)

		sfigure(2); clf;
		plotControl(p, a, u)
		% pause
	end % if p.plotLevel > 1

	if p.storeLevel
		U(:,p.nt+1-i) = u;
		V(:,p.nt+1-i) = v;
		% Ind(:,p.nt+1-i) = ind;
		% Vy(:,p.nt+1-i) = dv.dvdy;
	end

	v_prev = v;
	p.global_iter = p.global_iter + 1;

	if p.writeLevel > 1
		exportVTKSolution(p, a, v_prev, u, sprintf('solution_%06d.vtk',p.global_iter) );
	end % if p.writeLevel > 0

	if p.postprocessLevel == 2
		[dec, h] = postprocessSolution(p, c, a, v, 1);
		title(sprintf('Optimal control at t = %4d',p.t))
		saveas(h, sprintf([p.outputdir, '/', p.prefix, 'postprocessedSolution_%06d.png'], p.nt+1-i))
	end % if p.postprocessLevel == 2
end
