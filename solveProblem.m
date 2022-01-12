% This script solves the different HJB equations with post-processing
clear all, clc
% prob = 0; % Minimum Arrival Time 1d
% prob = 1; % Minimum Arrival Time 2d
% prob = 2; % Minimum Arrival Time 3d
% prob = 3; % Energy Storage 2d
prob = 4; % Energy Storage 3d
% prob = 5; % Worst of two asset put
% prob = 6; % Asian call option
% prob = 7; % American worst of two assets put option


% Start timer
tic

% Load helper functions such as positive, negative part
helperFunctions

p.method = 'FirstOrderOptimality';
switch prob
case 0
	p.dim = 1;
	p.mode = 'MinimumArrivalTime';
	p.outputdir = 'output1d';
	p.trivialControl = 0;
    
    p = example_MinimumArrivalTime_Problem(p);

case 1
	p.dim = 2;
	p.mode = 'MinimumArrivalTime';
	p.outputdir = 'output2d';
	p.trivialControl = 0;
    
    p = example_MinimumArrivalTime_Problem(p);

case 2
	p.dim = 3;
	p.mode = 'MinimumArrivalTime';
	p.outputdir = 'output3d';
	p.trivialControl = 0;

    p = example_MinimumArrivalTime_Problem(p);

case 3
	p.dim = 2;
	p.mode = 'EnergyStorage';
	p.outputdir = 'output2d';
	p.trivialControl = 0;

    p.method = 'SemiLagrangian';
    p = example_EnergyStorage_Problem(p);

case 4
	p.dim = 3;
	p.mode = 'EnergyStorage';
	p.outputdir = 'output3d';
	p.trivialControl = 0;

    % p.method = 'SemiLagrangian';
    p = example_EnergyStorage_Problem(p);

case 5
	p.dim = 2;
	p.mode = 'WorstOfTwoAssetCall';
	p.outputdir = 'tmp';
	p.trivialControl = 1;

    p = example_WorstAssetCall_Problem(p);

case 6
	p.dim = 2;
	p.mode = 'AsianCall';
	p.outputdir = 'tmp';
	p.trivialControl = 1;

    % TODO Implement semi lagrangian scheme for this case
    % p.method = 'SemiLagrangian';
    p = example_AsianCall_Problem(p);

case 7
	p.dim = 2;
	p.mode = 'AmericanWorstOfTwoAssetPut';
	p.outputdir = 'tmp';
	p.trivialControl = 0;

    p = example_AmericanWorstAssetPut_Problem(p);
end % switch prob

p.prefix = 'test';
% p.prefix = 'cv_final_35_'
% p.prefix = 'SemiLagrange_20_'

% Set some parameters
% p.method = 'SemiLagrangian'
p.parallel = 0;
p.refinement_level = 0;
% p.theta = 0.5;
p.theta = 1.0;
p.howard_maxiter = 50;
p.useUpwinding = 1;
p.convectionExplicit = 0;
p.damping = 0;

p.plotLevel  = 1;
p.writeLevel = 0;
p.storeLevel = 0;
p.postprocessLevel = 0;
p.errorLevelAbsolute = 1e-8;

% THIS IS THE POINT, WHERE COMPAREORIGANDPOSTCONTROL.M VARIES OVER DIFFERENT PARAMETER SETTINGS

% c = setupCoefficients(p);
[a, p] = setupStructure(p);

if p.postprocessLevel > 0
	p = initPostprocessing(p, a);
end

if strcmp(p.mode, 'EnergyStorage')
	p.Umin = p.alpha(a.XY);
	p.Umax = p.beta(a.XY);
    p.Uzero = zeros(p.ndofs, 1);
    p.Uzero(p.Umin>0) = p.Umin(p.Umin>0);
	
	if strcmp(p.method, 'SemiLagrangian')
		a = setupSemiLagrangian(p, a);
	end

end % if strcmp(p.mode, 'EnergyStorage')

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

toc
% Store workspace for further investigation
% save(['energy3d_', p.prefix, '.mat'])


if p.plotLevel > 0
	sfigure(1); clf;
	plotValue(p, a, v)

	sfigure(2); clf;
	plotControl(p, a, u)
	if p.dim == 3
		sfigure(3); clf;
		[X,Y,Z] = meshgrid(a.x, a.y, a.z);
		slice(X,Y,Z,p.vec2Mat(v),.5*(p.xmax-p.xmin),.5*(p.ymax-p.ymin),[.2 .8]*(p.zmax-p.zmin));

	end % if p.dim == 3
    if isfield(p, 'solution')
        sfigure(3); clf;
        vstar = p.solution(p.t, a.XY);
	    plotValue(p, a, vstar - v)
    end

end % if p.plotLevel > 0

if p.writeLevel > 0
	exportVTKSolution(p, a, v, u, 'solutionFinal.vtk');
end % if p.writeLevel > 0

if p.postprocessLevel == 1
	dec = postprocessSolution(p, c, a, v, 1);
end % if p.postprocessLevel == 1

% if p.postprocessLevel == 2
% 	if p.dim == 2 & strcmp(p.mode, 'EnergyStorage')
% 
% 		sfigure(3); clf, hold on
% 		% for i = 1:length(decBuy)
% 		for j = 1:20
% 			t = a.t(end -j + 1);
% 
% 			decBnd = decBuy{j};
% 			decBnd(:,3) = t;
% 			plot3(reshape(decBnd(:,1),2,[]), reshape(decBnd(:,2),2,[]), reshape(decBnd(:,3),2,[]), 'k-');
% 			decBnd = decSell{j};
% 			decBnd(:,3) = t;
% 			plot3(reshape(decBnd(:,1),2,[]), reshape(decBnd(:,2),2,[]), reshape(decBnd(:,3),2,[]), 'r-');
% 
% 		end % for i = 1:length(decBuy)
% 
% 	else
% 		fprintf('\n\nPostprocessing for 3d problem and minimum arrival time is not yet implemented\n\n')
% 	end % if p.dim == 2 & strcmp(p.mode, 'EnergyStorage')
% end % if p.postprocessLevel == 2
% 
