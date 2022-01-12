% This script solves the different HJB equations with post-processing
clear all, clc
set(0,'defaultAxesFontSize',24)
% set(0,'DefaultFigureVisible','off')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% 11 vs 15

% prob = 0; % 1d eikonal equation
% prob = 1; % 2d eikonal equation
% prob = 2; % 3d eikonal equation
% prob = 3; % 1d no quadr cost equation
% prob = 4; % 2d no quadr cost equation
% prob = 5; % 3d no quadr cost equation
% prob = 6; % 1d original
% prob = 7; % 2d original
% prob = 8; % 3d original
% prob = 9; % 2d Europ. WorstAssetPut with vanilla bc's
% prob = 10; % 2d Europ. WorstAssetCall with vanilla bc's
% prob = 11; % 2d Europ. WorstAssetPut with Diss bc's
% prob = 12; % 2d Europ. WorstAssetCall with Diss bc's
% prob = 13; % 2d Amer. WorstAssetPut with Diss bc's refined
% prob = 14; % 2d Amer. WorstAssetCall with Diss bc's
% prob = 15; % 2d Europ. WorstAssetPut with Diss bc's, refined at center
% prob = 16; % 2d Asian option with Howard
% prob = 17; % 2d Asian option with semi-Lagrange
% prob = 18; % 3d Energy Storage with Howard
prob = 19; % 3d Energy Storage with semi-Lagrange
% prob = 20; % 3d Energy Storage refined with Howard
% prob = 21; % 3d Energy Storage refined with semi-Lagrange

% Start timer
tic

% Load helper functions such as positive, negative part
helperFunctions

p.method = 'FirstOrderOptimality';
p.outputdir = 'output_MAT';
p.plotLevel  = 0;
p.postprocessLevel = 0;
p.trivialControl = 0;
p.refinement.level = 0;
p.solver = 'Howard';
problem_setup = @(p) example_MinimumArrivalTime_Problem(p);
nx0 = 1;

switch prob
case 0 % 1d eikonal equation
	p.dim = 1;

	p.trivialControl = 0;
    p.alpha = 0.0;
    p.beta = 0.0;
    p.sigmax = 0.0;
    p.xmax = +1.0;
    p.xmin = -1.0;
    Lset = 7:13;
    Nx = Lset;
    p.prefix = sprintf('experiment_MAT_Eikonal_dim_%d', p.dim);
    ntfunc = @(d, nx, l) d*nx;
    hfunc = @(dx, dt) dx + dt;
    % hfunc = @(dx,dt) sqrt(dx.^2 + dt)
    % ntfunc = @(d, nx, l) d*nx.^2;

    p.solution = @(x) 1-max(abs(x),[],2);   
    % Sol if xmin == 0
    % p.solution = @(x) min([1-x,x],[],2;

case 1 % 2d eikonal equation
	p.dim = 2;

	p.trivialControl = 0;
    p.alpha = 0.0;
    p.beta = 0.0;
    p.sigmax = 0.0;
    p.sigmay = 0.0;
    p.xmax = +1.0;
    p.xmin = -1.0;
    p.ymin = -1.0;
    Lset = 4:8;
    p.prefix = sprintf('experiment_MAT_Eikonal_dim_%d', p.dim);
    ntfunc = @(d, nx, l) d*nx;
    hfunc = @(dx,dt) dx + dt;
    p.solution = @(x) 1-max(abs(x),[],2);   

case 2 % 3d eikonal equation
	p.dim = 3;

	p.trivialControl = 0;
    p.alpha = 0.0;
    p.beta = 0.0;
    p.sigmax = 0.0;
    p.sigmay = 0.0;
    p.sigmaz = 0.0;
    p.xmax = +1.0;
    p.xmin = -1.0;
    p.ymin = -1.0;
    p.zmin = -1.0;
    Lset = 2:6;
    p.prefix = sprintf('experiment_MAT_Eikonal_dim_%d', p.dim);
    ntfunc = @(d, nx, l) d*nx;
    hfunc = @(dx,dt) dx + dt
    p.solution = @(x) 1-max(abs(x),[],2);   

case 3
	p.dim = 1;

	p.trivialControl = 0;
    p.sigmax = 0.5;
    p.xmax = +1.0;
    p.xmin = -1.0;

    p.alpha = 0.0;
    p.beta = 0.5;
    ntfunc = @(d, nx, l) 1*nx.^2;
    % ntfunc = @(d, nx, l) 2*nx;
    hfunc = @(dx,dt) sqrt(dx.^2 +dt)

    Lset = 2:8;
    p.prefix = sprintf('experiment_MAT_NoQuadrCost_dim_%d', p.dim);
    Nx = Lset;

case 4
	p.dim = 2;

	p.trivialControl = 0;
    p.sigmax = 0.5;
    p.sigmay = 0.5;
    p.xmax = +1.0;
    p.xmin = -1.0;
    p.ymax = +1.0;
    p.ymin = -1.0;

    p.alpha = 0.0;
    p.beta = 0.5;
    ntfunc = @(d, nx, l) 1*nx.^2;
    hfunc = @(dx,dt) sqrt(dx.^2 +dt);
    % ntfunc = @(d, nx, l) 2*nx;

    Lset = 1:7;
    p.prefix = sprintf('experiment_MAT_NoQuadrCost_dim_%d', p.dim);
    Nx = Lset;

case 5
	p.dim = 3;

	p.trivialControl = 0;
    p.sigmax = 0.5;
    p.sigmay = 0.5;
    p.sigmaz = 0.5;
    p.xmax = +1.0;
    p.xmin = -1.0;
    p.ymax = +1.0;
    p.ymin = -1.0;
    p.zmax = +1.0;
    p.zmin = -1.0;

    p.alpha = 0.0;
    p.beta = 0.5;
    ntfunc = @(d, nx, l) 1*nx.^2;
    hfunc = @(dx,dt) sqrt(dx.^2 +dt)
    % ntfunc = @(d, nx, l) 2*nx;

    Lset = 1:6;
    p.prefix = sprintf('experiment_MAT_NoQuadrCost_dim_%d', p.dim);
    Nx = Lset;

case 6
	p.dim = 1;

	p.trivialControl = 0;
    p.sigmax = 0.5;
    p.xmax = +1.0;
    p.xmin = -1.0;

    p.alpha = 0.5;
    p.beta = 0.5;
    ntfunc = @(d, nx, l) 1*nx.^2;
    % ntfunc = @(d, nx, l) 2*nx;
    hfunc = @(dx,dt) sqrt(dx.^2 +dt)

    Lset = 4:13;
    p.prefix = sprintf('experiment_MAT_Orig_dim_%d', p.dim);
    Nx = Lset;

case 7
	p.dim = 2;

	p.trivialControl = 0;
    p.sigmax = 0.5;
    p.sigmay = 0.5;
    p.xmax = +1.0;
    p.xmin = -1.0;
    p.ymax = +1.0;
    p.ymin = -1.0;

    p.alpha = 0.5;
    p.beta = 0.5;
    ntfunc = @(d, nx, l) 1*nx.^2;
    hfunc = @(dx,dt) sqrt(dx.^2 +dt)
    % ntfunc = @(d, nx, l) 2*nx;

    Lset = 1:7;
    p.prefix = sprintf('experiment_MAT_Orig_dim_%d', p.dim);
    Nx = Lset;

case 8
	p.dim = 3;

	p.trivialControl = 0;
    p.sigmax = 0.5;
    p.sigmay = 0.5;
    p.sigmaz = 0.5;
    p.xmax = +1.0;
    p.xmin = -1.0;
    p.ymax = +1.0;
    p.ymin = -1.0;
    p.zmax = +1.0;
    p.zmin = -1.0;

    p.alpha = 0.5;
    p.beta = 0.5;
    ntfunc = @(d, nx, l) 1*nx.^2;
    hfunc = @(dx,dt) sqrt(dx.^2 +dt);
    % ntfunc = @(d, nx, l) 2*nx;

    Lset = 1:6;
    p.prefix = sprintf('experiment_MAT_Orig_dim_%d', p.dim);
    Nx = Lset;

case 9
	p.dim = 2;
    nx0 = 10;

	p.dim = 2;
    p.outputdir = 'output_WorstAssetOption';
	p.trivialControl = 1;
    p.payoff = 'put';
    p.type = 'European';
    p.bctype = 1;
    p.prefix = sprintf('experiment_WorstAssetPut_dim_%d', p.dim);
    problem_setup = @(p) example_WorstAssetOption_Problem(p);
    ntfunc = @(d, nx, l) round(.1*nx.^2)+1;
    hfunc = @(dx,dt) sqrt(dx.^2 +dt);
    % ntfunc = @(d, nx, l) 2*nx;

    Lset = 0:5;
    Nx = Lset;

case 10
    nx0 = 10;

	p.dim = 2;
    p.outputdir = 'output_WorstAssetOption';
	p.trivialControl = 1;
    p.payoff = 'call';
    p.type = 'European';
    p.bctype = 1;
    p.prefix = sprintf('experiment_WorstAssetCall_dim_%d', p.dim);
    problem_setup = @(p) example_WorstAssetOption_Problem(p);
    ntfunc = @(d, nx, l) round(0.1*nx.^2);
    hfunc = @(dx,dt) sqrt(dx.^2 +dt);
    % ntfunc = @(d, nx, l) 2*nx;

    Lset = 2:6;
    Nx = Lset;

case 11

	p.dim = 2;
    p.outputdir = 'output_WorstAssetOption';
	p.trivialControl = 1;
    p.payoff = 'put';
    p.type = 'European';
    p.bctype = 2;
    p.prefix = sprintf('experiment_WorstAssetPut_bc2_dim_%d', p.dim);
    problem_setup = @(p) example_WorstAssetOption_Problem(p);
    % Works also for explicit computations
    ntfunc = @(d, nx, l) round(000.5*nx.^2)+1;
    hfunc = @(dx,dt) sqrt(dx.^2 +dt);
    % Only possible for theta = 1
    ntfunc = @(d, nx, l) l*nx;
    hfunc = @(dx,dt) dx +dt;
    nx0 = 7;
    p.xmax = 140;
    p.ymax = 140;
    p.xlabel = '$x$';
    p.ylabel = '$y$';

    % ntfunc = @(d, nx, l) 2*nx;
    Lset = 1:6;

    % Lset = 2:6;
    Nx = Lset;

case 12
	p.dim = 2;
    nx0 = 10;

	p.dim = 2;
    p.outputdir = 'output_WorstAssetOption';
	p.trivialControl = 1;
    p.payoff = 'call';
    p.type = 'European';
    p.bctype = 2;
    p.prefix = sprintf('experiment_WorstAssetCall_bc2_dim_%d', p.dim);
    problem_setup = @(p) example_WorstAssetOption_Problem(p);
    ntfunc = @(d, nx, l) round(0.1*nx.^2);
    hfunc = @(dx,dt) sqrt(dx.^2 +dt);
    % ntfunc = @(d, nx, l) 2*nx;
    p.xlabel = '$x$';
    p.ylabel = '$y$';

    Lset = 2:6;
    Nx = Lset;

case 13
	p.dim = 2;
    nx0 = 10;

	p.dim = 2;
    p.outputdir = 'output_WorstAssetOption';
    p.payoff = 'put';
    p.type = 'American';
    p.bctype = 2;
    p.prefix = sprintf('experiment_AmericanWorstAssetPut_bc2_dim_%d', p.dim);
    problem_setup = @(p) example_WorstAssetOption_Problem(p);
    p.xlabel = '$x$';
    p.ylabel = '$y$';
    
    % ntfunc = @(d, nx, l) round(000.5*nx.^2)+1;
    % hfunc = @(dx,dt) sqrt(dx.^2 +dt);

    % Only possible for theta = 1
    ntfunc = @(d, nx, l) 5*nx+1;
    hfunc = @(dx,dt) dx +dt;
    p.refinement.level = 1;
    p.refinement.mode = 'abs';
    p.refinement.xlow = [20,30];
    p.refinement.xup = [60,50];
    p.refinement.ylow = [20,30];
    p.refinement.yup = [60,50];
    p.prefix = sprintf('experiment_AmericanWorstAssetPut_bc2_refCenter_dim_%d', p.dim);

    nx0 = 14;
    p.xmax = 140;
    p.ymax = 140;

    Lset = 1:5;
    % Lset = 2:6;
    Nx = Lset;

case 14
	p.dim = 2;
    p.outputdir = 'output_WorstAssetOption';
    p.payoff = 'call';
    p.type = 'American';
    p.bctype = 2;
    p.prefix = sprintf('experiment_AmericanWorstAssetCall_bc2_dim_%d', p.dim);
    problem_setup = @(p) example_WorstAssetOption_Problem(p);
    ntfunc = @(d, nx, l) round(0.1*nx.^2);
    hfunc = @(dx,dt) sqrt(dx.^2 +dt);
    % ntfunc = @(d, nx, l) 2*nx;
    nx0 = 14;
    p.xmax = 140;
    p.ymax = 140;
    p.xlabel = '$x$';
    p.ylabel = '$y$';

    % Lset = 1:5;
    Lset = 2:6;
    Nx = Lset;

case 15
	p.dim = 2;
    p.outputdir = 'output_WorstAssetOption';
	p.trivialControl = 1;
    p.payoff = 'put';
    p.type = 'European';
    p.bctype = 2;
    p.prefix = sprintf('experiment_WorstAssetPut_bc2_refCenter_dim_%d', p.dim);
    problem_setup = @(p) example_WorstAssetOption_Problem(p);
    % Works also for explicit computations
    ntfunc = @(d, nx, l) round(0.5*nx.^2)+1;
    hfunc = @(dx,dt) sqrt(dx.^2 +dt);
    % Only possible for theta = 1
    ntfunc = @(d, nx, l) l*nx;
    hfunc = @(dx,dt) dx +dt;
    nx0 = 7;
    p.xmax = 140;
    p.ymax = 140;
    p.refinement.level = 1;
    p.refinement.mode = 'abs';
    p.refinement.xlow = [20,30];
    p.refinement.xup = [60,50];
    p.refinement.ylow = [20,30];
    p.refinement.yup = [60,50];
    p.xlabel = '$x$';
    p.ylabel = '$y$';

    % ntfunc = @(d, nx, l) 2*nx;
    Lset = 1:6;

    % Lset = 2:6;
    Nx = Lset;

case 16
    % Coarse 1.7611, finer 2.0564
	p.dim = 2;
    p.outputdir = 'output_AsianOption';
	p.trivialControl = 1;
    % p.payoff = 'put';
    % p.type = 'European';
    p.bctype = 2;
    p.prefix = sprintf('experiment_AsianCall');
    problem_setup = @(p) example_AsianCall_Problem(p);
    % Works also for explicit computations
    % Only possible for theta = 1
    ntfunc = @(d, nx, l) nx;
    hfunc = @(dx,dt) dx + dt;
    nx0 = 4;
    p.xlabel = '$x$';
    p.ylabel = '$y$';
    p.refinement.level = 1;
    p.refinement.mode = 'abs';
    p.refinement.xlow = [90 80];
    p.refinement.xup = [110 120];
    p.refinement.ylow = [90 80];
    p.refinement.yup = [110 120];

    % ntfunc = @(d, nx, l) 2*nx;
    Lset = 1:7;

    % Lset = 2:6;
    Nx = Lset;

case 17

	p.dim = 2;
    p.outputdir = 'output_AsianOption';
	p.trivialControl = 1;
    p.solver = 'SL';
    p.prefix = sprintf('experiment_AsianCall_SL');
    problem_setup = @(p) example_AsianCall_Problem(p);
    % Works also for explicit computations
    % Only possible for theta = 1
    ntfunc = @(d, nx, l) nx;
    hfunc = @(dx,dt) dx + dt;
    p.xlabel = '$x$';
    p.ylabel = '$y$';
    nx0 = 4;
    p.refinement.level = 1;
    p.refinement.mode = 'abs';
    p.refinement.xlow = [90 80];
    p.refinement.xup = [110 120];
    p.refinement.ylow = [90 80];
    p.refinement.yup = [110 120];

    % ntfunc = @(d, nx, l) 2*nx;
    Lset = 1:8;

    % Lset = 2:6;
    Nx = Lset;

case 18

	p.dim = 3;
    p.outputdir = 'output_Energy';
    set(0,'defaultAxesFontSize',18)
    p.solver = 'Howard';
    p.prefix = sprintf('experiment_Energy_Howard_uniform_1-5_Theta_00');
    problem_setup = @(p) example_EnergyStorage_Problem(p);
    p.Tmin = 0;
    ntfunc = @(d, nx, l) l*(365-p.Tmin);
    hfunc = @(dx,dt) dx + dt;
    nx0 = 2;
    p.xlabel = '$p$';
    p.ylabel = '$q$';

    Lset = 1:5;
    p.postprocessLevel = 0;
    Nx = Lset;

case 19

	p.dim = 3;
    p.outputdir = 'output_Energy';
    set(0,'defaultAxesFontSize',18)
    p.solver = 'SL';
    p.prefix = sprintf('experiment_Energy_SL_uniform_1-5');
    problem_setup = @(p) example_EnergyStorage_Problem(p);
    p.Tmin = 0;
    ntfunc = @(d, nx, l) l*(365-p.Tmin);
    hfunc = @(dx,dt) dx + dt;
    nx0 = 2;
    Lset = 1:5;
    p.postprocessLevel = 0;
    Nx = Lset;
    p.xlabel = '$p$';
    p.ylabel = '$q$';

case 20

	p.dim = 3;
    p.outputdir = 'output_Energy';
    set(0,'defaultAxesFontSize',18)
    p.solver = 'Howard';
    problem_setup = @(p) example_EnergyStorage_Problem(p);
    p.Tmin = 0;
    ntfunc = @(d, nx, l) l*(365-p.Tmin);
    hfunc = @(dx,dt) dx + dt;
    nx0 = 2;
    p.prefix = sprintf('experiment_Energy_Howard_ref');
    p.refinement.level = 1;
    p.refinement.mode = 'abs';
    p.refinement.xlow = [30 20];
    p.refinement.xup = [50 60];
    p.refinement.ylow = [80000];
    p.refinement.yup = [140000];
    p.refinement.zlow = [7000];
    p.refinement.zup = [9000];
    p.xlabel = '$p$';
    p.ylabel = '$q$';

    Lset = 1:4;
    p.postprocessLevel = 0;
    Nx = Lset;

case 21

	p.dim = 3;
    p.outputdir = 'output_Energy';
    set(0,'defaultAxesFontSize',18)
    p.solver = 'SL';
    problem_setup = @(p) example_EnergyStorage_Problem(p);
    p.Tmin = 0;
    ntfunc = @(d, nx, l) l*(365-p.Tmin);
    hfunc = @(dx, dt) dx + dt;
    nx0 = 2;
    p.xlabel = '$p$';
    p.ylabel = '$q$';
    p.prefix = sprintf('experiment_Energy_SL_ref');
    p.refinement.level = 1;
    p.refinement.mode = 'abs';
    p.refinement.xlow = [30 20];
    p.refinement.xup = [50 60];
    p.refinement.ylow = [80000];
    p.refinement.yup = [140000];
    p.refinement.zlow = [7000];
    p.refinement.zup = [9000];

    Lset = 1:4;
    p.postprocessLevel = 0;
    Nx = Lset;

end

Nx = Lset;
Ndofs = (nx0*2.^Nx+1).^p.dim;

% Thetaset = [0.0, 0.5, 1.0];
% fheader={'Level', 'N_dofs', 'N_t', 'Delta_x', 'Delta_t', 'h', 't0_0abs', 't0_0p', 't0_5abs', 't0_5p', 't1_0abs', 't1_0p', 'qoi_t0_0abs', 'qoi_t0_0p', 'qoi_t0_0val', 'qoi_t0_5abs', 'qoi_t0_5p', 'qoi_t0_5val', 'qoi_t1_0abs', 'qoi_t1_0p', 'qoi_t1_0val'};
Thetaset = [1.0];
fheader={'Level', 'N_dofs', 'N_t', 'Delta_x', 'Delta_t', 'h', 'time', 't1_0abs', 't1_0p', 'qoi_t1_0abs', 'qoi_t1_0p', 'qoi_t1_0val'};
TL = length(Thetaset);
off1 = 7;
off2 = off1+2*TL;

Out = zeros(length(Lset),length(Thetaset)*2+2);
Out(:,1) = Lset;
k = 1;
for theta = Thetaset
    lidx = 1;
    for l = Lset

        if l == 3 & theta == 1
            p.writeValueFunction = 1;
            p.writeControlFunction = 1;
        else
            p.writeValueFunction = 0;
            p.writeControlFunction = 0;
        end

        p.theta = theta;
        p.nx = nx0*2^l;
        if p.dim > 1
            p.ny = nx0*2^l;
        end
        if p.dim > 2
            p.nz = nx0*2^l;
        end
        % p.nt = nt0*2^l;
        p.nt = ntfunc(p.dim, p.nx, l);
        p = problem_setup(p);
    
        
        % Set some parameters
        p.parallel = 0;
        p.howard_maxiter = 50;
        p.useUpwinding = 1;
        p.convectionExplicit = 0;
        p.damping = 0;
        
        p.writeLevel = 0;
        p.storeLevel = 0;
        p.errorLevelAbsolute = 1e-10;
        
        % c = setupCoefficients(p);
        [a, p] = setupStructure(p);

        Out(lidx,2) = p.ndofs;
        Out(lidx,3) = p.nt;
        
        a.norm1 = @(v) norm_Linfty(v);
        a.norm2 = @(v) norm_L2(a.Hv, v);
        
        if p.plotLevel > 3
        	sfigure(3); clf
        	plotDiffusion(p, a);
        	sfigure(4); clf
        	plotConvection(p, a, 0, zeros(size(a.X)));
            keyboard
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

        if p.postprocessLevel > 0
            p = initPostprocessing(p, a);
        end
        
        % Start temporal loop
        for i = 1:p.nt
        
        	% Set current time
        	p.t = a.t(p.nt+1-i);
            switch p.solver
                case 'Howard'
                    [v, u] = solveHowardIteration(p, a, v_prev, u_prev);
                case 'SL'
                    if i == 1
                        a = setupSemiLagrangian(p, a);
                    end
                    [v, u] = solveSemiLagrangeIteration(p, a, v_prev, u_prev);
            end

        
        	if p.plotLevel > 1
        		sfigure(1); clf;
        		plotValue(p, a, v)
                if ~p.trivialControl
            		sfigure(2); clf;
            		plotControl(p, a, u)
                end
                pause
        	end % if p.plotLevel > 1
        
        	if p.storeLevel
        		U(:,p.nt+1-i) = u;
        		V(:,p.nt+1-i) = v;
        		% Ind(:,p.nt+1-i) = ind;
        		% Vy(:,p.nt+1-i) = dv.dvdy;
            end

            if isfield(p, 'postprocessSolution')
                if p.postprocessSolution
                    dec = postprocessSolution(p, a, v, 1);
                end
            end
        
        	v_prev = v;
        	p.global_iter = p.global_iter + 1;
        
        	if p.writeLevel > 1
        		exportVTKSolution(p, a, v_prev, u, sprintf('solution_%06d.vtk',p.global_iter) );
        	end % if p.writeLevel > 0
        
        	if p.postprocessLevel == 2
        		[dec, h] = postprocessSolution(p, a, v, 1);
        		title(sprintf('Optimal control at t = %5.2f',p.t))
                P = 20;
                set(gcf, 'Units', 'inches');
                set(gcf, 'PaperPosition', [0 0 P P]); %Position plot at left hand corner with width 5 and height 5.
                set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
                fname = sprintf([p.outputdir, '/', p.prefix, '_postprocessedSolution_%06d.png'], p.nt+1-i);
                print(gcf,fname,'-dpng','-r300'); 
        	end % if p.postprocessLevel == 2
        end
        
        trun = toc;
        fprintf('Time for current run (s): %6.2f', trun)
        Out(lidx, 7) = trun;
        switch p.dim
        case 1
            Vq{l+1} = @(x) interp1(a.x,v,x);
        case 2
            ax = p.vec2Mat(a.X);
            ay = p.vec2Mat(a.Y);
            av = p.vec2Mat(v);
            Vq{l+1} = @(x) interp2(ax, ay, av, x(:,1), x(:,2));
        case 3
            ax = p.vec2Mat(a.X);
            ay = p.vec2Mat(a.Y);
            az = p.vec2Mat(a.Z);
            av = p.vec2Mat(v);
            Vq{l+1} = @(x) interp3(ax, ay, az, av, x(:,1), x(:,2), x(:,3));
        end

        
        fprintf('Value in qoi: %6.4f\n', Vq{l+1}(p.qoi))
        
        if p.plotLevel > 0
        	sfigure(1); clf;
            plotValue(p, a, v)
       
            if ~p.trivialControl
                sfigure(2); clf;
            	plotControl(p, a, u)
            	if p.dim == 3
            		sfigure(3); clf;
            		[X, Y, Z] = meshgrid(a.x, a.y, a.z);
            		slice(X, Y, Z, p.vec2Mat(v), .5*(p.xmax-p.xmin),.5*(p.ymax-p.ymin),[.2 .8]*(p.zmax-p.zmin));
            
            	end % if p.dim == 3
            end
            if isfield(p, 'solution')
                sfigure(3); clf;
                % vstar = p.solution(p.t, a.XY);
                vstar = p.solution(a.XY);
                plotValue(p, a, vstar - v, [], 'comparison');
        	    plotValue(p, a, abs(vstar - v), [], 'comparison_abs');
                plotValue(p, a, abs(vstar - v)./vstar, [], 'relError');
            end
            keyboard
        
        end % if p.plotLevel > 0
        
        if p.writeLevel > 0
        	exportVTKSolution(p, a, v, u, 'solutionFinal.vtk');
        end % if p.writeLevel > 0
    
        switch p.dim
        case 1
            Deltax = max(a.dx);
        case 2
            Deltax = max([a.dx,a.dy]);
        case 3
            Deltax = max([a.dx,a.dy,a.dz]);
        end
        Deltat = max(diff(a.t));
        Out(lidx,4) = Deltax;
        Out(lidx,5) = Deltat;
        Out(lidx,6) = sqrt(Deltax^2 + Deltat);
        lidx = lidx + 1;
    end % for l=Lset
    
    i = 1;
    for l=Lset
        if isfield(p, 'solution')
            vdiff = Vq{l+1}(a.XY)-p.solution(a.XY);
            v_qoi = Vq{l+1}(p.qoi);
            vdiff_qoi = v_qoi - p.solution(p.qoi);
        else
            vdiff = v - Vq{l+1}(a.XY);
            v_qoi = Vq{l+1}(p.qoi);
            vdiff_qoi = v_qoi - Vq{Lset(end)+1}(p.qoi);
        end
        vdiff_Linfty_norm(i) = norm_Linfty(vdiff);
        v_qoi_val(i) = v_qoi;
        vdiff_qoi_Linfty_norm(i) = norm_Linfty(vdiff_qoi);
        fprintf('Error level %d: %6.2e\n', l, vdiff_Linfty_norm(i));
        i = i+1;
    end
    
    sfigure(5); clf
    %Nx = 1./(nx0*2.^((1:(L-1))-1));
    Vx = vdiff_Linfty_norm;
    plot(Nx, log(Vx)/log(2),'k+-');

    fprintf('p_est w.r.t total dofs\n')
    Tdofs = (Out(:,2).*Out(:,3))';
    p_est = -(p.dim+1)*log(Vx(2:end)./Vx(1:end-1))/log(Tdofs(2:end)/Tdofs(1:end-1))

    fprintf('p_est w.r.t spatial dofs\n')
    p_est = -p.dim*log(Vx(2:end)./Vx(1:end-1))/log(Ndofs(2:end)/Ndofs(1:end-1))
    % hval = Ndofs
    % hval = sqrt()Ndofs
    % p_est = -p.dim*log(Vx(2:end)./Vx(1:end-1))/log(hval(2:end)/hval(1:end-1))
    hval = hfunc(Out(:,4)',Out(:,5)')
    fprintf('p_est w.r.t hval\n')
    p_est = log(Vx(2:end)./Vx(1:end-1))/log(hval(2:end)/hval(1:end-1))
    Out(:,off1+1+2*(k-1)) = Vx;
    Out(2:end,off1+2+2*(k-1))=p_est;

    % Estimate values at quantity of interest
    Vx = vdiff_qoi_Linfty_norm;
    Out(:,off2+1+3*(k-1)) = Vx;
    p_est = log(Vx(2:end)./Vx(1:end-1))/log(hval(2:end)/hval(1:end-1))
    Out(2:end,off2+2+3*(k-1)) = p_est;
    Out(:,off2+3+3*(k-1)) = v_qoi_val;
    k = k+1;


    % Write csv file with function values on finest mesh
    % if p.dim == 2
    %     if a.nx == 32
    %         fname = strcat(p.outputdir,'/', p.prefix, sprintf('_value_theta_%f.csv',p.theta));
    %         fheader={'x', 'y', 'z'};
    %         fid = fopen(fname, 'w');
    %         fprintf(fid, '%s\n', strjoin(fheader,','));
    %         D = [a.X,a.Y,v];
    %         fclose(fid)
    %         dlmwrite(fname, D(1:(a.nx+1),:), '-append');
    %         fid = fopen(fname, 'a');
    %         fprintf(fid, '\n');
    %         dlmwrite(fname, D((a.nx+2):end,:), '-append');
    %         fclose(fid)

    %         % Plot control
    %         fname = strcat(p.outputdir,'/', p.prefix, sprintf('_control_theta_%f.csv',p.theta));
    %         fheader={'x', 'y', 'u1', 'u2'};
    %         fid = fopen(fname, 'w');
    %         fprintf(fid, '%s\n', strjoin(fheader,','));
    %         D = [a.X,a.Y,u];
    %         fclose(fid)
    %         dlmwrite(fname, D(1:(a.nx+1),:), '-append');
    %         fid = fopen(fname, 'a');
    %         fprintf(fid, '\n');
    %         dlmwrite(fname, D((a.nx+2):end,:), '-append');
    %         fclose(fid)
    %     end
    % end
end


% Linear regression to estimate order of convergence
Hx = Out(:,4);
if isfield(p,'solution')
    kl = length(Hx);
else
    kl = length(Hx)-2;
end

for k=1:TL
    Hxlog=log(hfunc(Out(:,4),Out(:,5)))';

    AA = [sum(Hxlog(1:kl).^2), sum(Hxlog(1:kl));sum(Hxlog(1:kl)), length(Hxlog(1:kl))];
    
    % Convergence order linfty norm by linear interpolation
    AAr = [sum(Hxlog(1:kl).*log(Out(1:kl,off1+1+2*(k-1)))'); sum(log(Out(1:kl,off1+1+2*(k-1))))];
    tmp = AA\AAr
    tp = tmp(1);
    Out(length(Hx)+1,off1+2+2*(k-1)) = tp;

    % Convergence order qoi by linear interpolation
    AAr = [sum(Hxlog(1:kl).*log(Out(1:kl,off2+1+3*(k-1)))'); sum(log(Out(1:kl,off2+1+3*(k-1))))];
    tmp = AA\AAr;
    tp = tmp(1);
    Out(length(Hx)+1,off2+2+3*(k-1)) = tp;
end

% Write csv with convergence results
fname = strcat(p.outputdir,'/', p.prefix, '.csv');
fid = fopen(fname, 'w');
fprintf(fid, '%s\n', strjoin(fheader,','));
dlmwrite(fname, Out, '-append','precision',10);

% Replace zeros by empty spaces
system(sprintf("sed -i  's/^0,/,/g' %s", fname));
system(sprintf("sed -i  's/,0$/,/g' %s", fname));
system(sprintf("sed -i  's/,0,/,,/g' %s", fname));
system(sprintf("sed -i  's/,0,/,,/g' %s", fname));
	
