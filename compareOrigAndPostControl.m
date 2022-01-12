% This script is derived from solveProblem.m and loops over different
% parameter settings to generate a table used for the illustration of 
% the postprocessing method in terms of the gas storage facility problem
clear all
%         np   nq   nc   nt  
PARAMS = [
					10   10   10  365;
					20   20   20  365;
					40   40   40  365;
					80   80   80  365;
					% 160   160   160  365
					% 10   10   10  730;
					% 20   20   20  730;
					% 40   40   40  730;
					% 80   80   80  730;
					% 120   60   60  365;
					% 10   10   10  730;
					% 20   20   20  730;
					% 40   40   40  730;
					% 80   80   80  730;
					% 60   60   60  365;
					% 100   100   100  365
					% 20   20   20  730;
					];

p.dim = 3;
p.mode = 'EnergyStorage';
p.outputdir = 'outputForPostprocessingPaper';
p.method = 'SemiLagrangian';
p.refinement_level = 1;
p.theta = 1.0;
p.howard_maxiter = 50;
p.useUpwinding = 1;
p.convectionExplicit = 0;
p.damping = 1.0;
p.plotLevel  = 0;
p.writeLevel = 0;
p.storeLevel = 1;
p.postprocessLevel = 0;
p.errorLevelAbsolute = 1e-8;

% Read parameters
p = setupParameters(p);

% Loop over parameters in PARAMS
for pi = 1:size(PARAMS,1)
	p.nx = PARAMS(pi,1);
	p.ny = PARAMS(pi,2);
	p.nz = PARAMS(pi,3);
	p.nt = PARAMS(pi,4);
	p.Tmax = p.nt

	% ----------------------------------------------------------------------------------------------
	
	c = setupCoefficients(p);
	[a, p] = setupStructure(p);
	if p.postprocessLevel > 0
		p = initPostprocessing(p, a);
	end
	
	if strcmp(p.mode, 'EnergyStorage')
		p.Umin = c.alpha(a.XY);
		p.Umax = c.beta(a.XY);
		
		if strcmp(p.method, 'CompareValues') && p.dim ==3
			p = setupCompareValues(p, a, c);
		end % if strcmp(p.method, 'CompareValues')
	
		if strcmp(p.method, 'SemiLagrangian')
			a = setupSemiLagrangian(p, a, c);
		end
	
	end % if strcmp(p.mode, 'EnergyStorage')
	
	
	a.norm1 = @(v) norm_Linfty(v);
	a.norm2 = @(v) norm_L2(a.Hv, v);
	
	if p.plotLevel > 3
		sfigure(3); clf
		plotDiffusion(p, c, a);
		sfigure(4); clf
		plotConvection(p, c, a, 0, zeros(size(a.X)));
		pause
	end
	
	% Setup second order matrix
	a.A2 = setupSecondOrderMatrix(p, c, a);
	
	% Setup routines to set up first order matrix/matrices
	[a.A1func] = setupFirstOrderMatrix(p, c, a);
	
	% Solve 
	% -v_t + sup_{u \in U} {-A(x): (nabla^2 v)(x) + B(x)' * \grad(v)(x) + C(x) - F(x)} = 0
	
	m = p.ndofs;
	n = p.ndofs;
	Id = speye(m);
	
	v_prev = c.finalTimeVal(a.XY);
	u_prev = zeros(size(v_prev));
	
	a.Ffunc = @(u) c.f(a.XY,u);
	
	
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
	
		% Print head of convergence table on each 10th iteration
		if mod(p.global_iter, 10) == 0
			fprintf('\n   Iter    |   Time   | ||w_n-w_nm1||_inf   |   ||w_n - w_nm1||_2  |   ||u_n - u_nm1||_2\n')
		end
		fprintf('---------------------------------------------------------------------------\n')
	
	
		if strcmp(p.method, 'SemiLagrangian')
			[v, u] = solveSemiLagrangeIteration(p, c, a, v_prev, u_prev);
		else
			[v, u] = solveHowardIteration(p, c, a, v_prev, u_prev);
		end
	
		if p.plotLevel > 1
			sfigure(1); clf;
			plotValue(p, c, a, v)
	
			sfigure(2); clf;
			plotControl(p, c, a, u)
			pause
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
	
	% Store workspace for further investigation
	% save(['energy3d_', p.prefix, '.mat'])
	
	
	if p.plotLevel > 0
		sfigure(1); clf;
		plotValue(p, c, a, v)
	
		sfigure(2); clf;
		plotControl(p, c, a, u)
		if p.dim == 3
			sfigure(3); clf;
			[X,Y,Z] = meshgrid(a.x, a.y, a.z);
			slice(X,Y,Z,p.vec2Mat(v),.5*(p.xmax-p.xmin),.5*(p.ymax-p.ymin),[.2 .8]*(p.zmax-p.zmin));
	
		end % if p.dim == 3
	
	end % if p.plotLevel > 0
	
	if p.writeLevel > 0
		exportVTKSolution(p, a, v, u, 'solutionFinal.vtk');
	end % if p.writeLevel > 0
	
	if p.postprocessLevel == 1
		dec = postprocessSolution(p, c, a, v, 1);
	end % if p.postprocessLevel == 1
	
	% ----------------------------------------------------------------------------------------------

	simulateValueVectorized;
	% keyboard
	% Write values-of-interest to PARAMS
	[X,Y,Z] = meshgrid(a.x, a.y, a.z);
	PARAMS(pi,5) = p.ndofs;
	PARAMS(pi,10) = interp3(X, Y, Z, p.vec2Mat(V(:,1)), x0(1), x0(2), x0(3));
	PARAMS(pi,11) = mean(F_Orig);
	PARAMS(pi,12) = mean(F_Post);
	PARAMS(pi,13) = mean(imp);
	PARAMS(pi,14) = std(imp);
	PARAMS(pi,15) = median(imp);
	drawnow
end

csvwrite([p.outputdir, '/compareOrigAndPostControl_w_corr_w_opportunity_w_cost_365_pmax80_refine.csv'], PARAMS);

