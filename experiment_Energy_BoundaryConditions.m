% This script solves the different HJB equations with post-processing
clear all, clc
% set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesFontSize',18)
% set(0,'DefaultFigureVisible','off')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


% prob = 18; % 3d Energy Storage with Howard
prob = 19; % 3d Energy Storage with semi-Lagrange
% prob = 20; % 3d Energy Storage refined with Howard
% prob = 21; % 3d Energy Storage refined with semi-Lagrange

% Start timer
tic

% Load helper functions such as positive, negative part
helperFunctions

p.method = 'FirstOrderOptimality';
p.outputdir = '/LOCAL/blja/EnergyDiss';
p.plotLevel  = 0;
p.postprocessLevel = 2;
p.trivialControl = 0;
p.refinement.level = 0;

Nxset = 32;
p.nt = 365*20;

% BC 3 nochmal mit diesen Einstellungen
% Nxset = 32;
% p.nt = 365*4;

p.Tmin = 0;
hfunc = @(dx,dt) dx + dt;
p.xlabel = '$p$';
p.ylabel = '$q$';
problem_setup = @(p) example_EnergyStorage_Problem(p);
p.dim = 3;
Lset = 0;
Nx = Lset;
p.bctype = 1;
savePostprocessSol = 1;

switch prob
case 18

    p.solver = 'Howard';
    p.prefix = sprintf('experiment_Energy_Howard_Expl_bc_%d',p.bctype);
    p.theta = 0.0;

case 19

    p.solver = 'SL';
    p.prefix = sprintf('experiment_Energy_SL_bc_%d',p.bctype);
    p.theta = 1.0;

case 20

    p.solver = 'Howard';
    p.prefix = sprintf('experiment_Energy_Howard_ref_Expl_bc_%d',p.bctype);
    p.refinement.level = 1;
    p.refinement.mode = 'abs';
    p.refinement.xlow = [30 20];
    p.refinement.xup = [50 60];
    p.refinement.ylow = [80000];
    p.refinement.yup = [140000];
    p.refinement.zlow = [7000];
    p.refinement.zup = [9000];
    p.theta = 0.0;

case 21

    p.solver = 'SL';
    p.prefix = sprintf('experiment_Energy_SL_ref_bc_%d',p.bctype);
    p.refinement.level = 1;
    p.refinement.mode = 'abs';
    p.refinement.xlow = [30 20];
    p.refinement.xup = [50 60];
    p.refinement.ylow = [80000];
    p.refinement.yup = [140000];
    p.refinement.zlow = [7000];
    p.refinement.zup = [9000];
    p.theta = 1.0;

end

k = 1;
for nx0 = Nxset
    Ndofs = nx0*ones(size(Lset));
    lidx = 1;
    for l = Lset

        p.writeValueFunction = 0;
        p.writeControlFunction = 0;

        p.nx = Ndofs(lidx);
        if p.dim > 1
            p.ny = Ndofs(lidx);
        end
        if p.dim > 2
            p.nz = Ndofs(lidx);
        end
        % p.nt = nt0*2^l;
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
                % pause
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
        
            if p.postprocessLevel == 2 & mod(p.t,1) == 0 | mod(p.t,1) == 0.5
            % if p.postprocessLevel == 2 & any(p.t==[364,363,362,361,360,330,300,270,240,210,180,150,120,90,60,30,0])
        		[dec, h] = postprocessSolution(p, a, v, 1);
        		title(sprintf('Optimal control at t = %5.2f',p.t));
                % pause
                % keyboard
                P = 20;
                set(gcf, 'Units', 'inches');
                set(gcf, 'PaperPosition', [0 0 P P]); %Position plot at left hand corner with width 5 and height 5.
                set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
                fname = sprintf([p.outputdir, '/', p.prefix, '_postprocessedSolution_t%06d.png'], p.t*100);
                if savePostprocessSol
                    print(gcf,fname,'-dpng','-r300'); 
                else
                    pause

                end
                
        	end % if p.postprocessLevel == 2
        end
        
        trun = toc;
        fprintf('Time for current run (s): %6.2f', trun)
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
        lidx = lidx + 1;
    end % for l=Lset
    
    k = k+1;


end
