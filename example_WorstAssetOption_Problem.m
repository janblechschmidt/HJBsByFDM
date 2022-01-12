function p = example_WorstAssetOption_Problem(p)
    
    n = @(x) size(x,1);
    EPS = 1e-12;

    if ~isfield(p, 'nt')
        p.nt = 320;
    end
    p.Tmin = 0.5;
    p.Tmax = 1;
    if ~isfield(p, 'nx')
        p.nx = 200;
    end
    p.xmin = 0;
    if ~isfield(p, 'xmax')
        p.xmax = 100;
    end
    if ~isfield(p, 'ny')
        p.ny = 200;
    end
    p.ymin = 0;
    if ~isfield(p, 'ymax')
        p.ymax = 100;
    end

    % tau = 0.5;
    H = 40;
    V = 40;
    F = 40;
    R = 0.05;
    sigmaH = 0.3;
    sigmaV = 0.3;
    rho_VH = 0.5;
    p.qoi = 40*ones(1,p.dim);
    p.control_dim = 2;

    p.potential    = @(x) R * ones(n(x),1);
    p.f            = @(x,u) zeros(n(x),1);

    p.DirichletBC = 1;

    % Put
    if ~isfield(p, 'payoff')
        p.payoff = 'put'
    end

    if ~isfield(p, 'type')
        p.type = 'European'
    end

    if ~isfield(p, 'bctype')
        p.bctype = 1;
    end

    switch p.payoff
    case 'put'
        p.finalTimeVal = @(x) max(F-min(x,[],2),0);
        if strcmp(p.type, 'European')
            p.boundaryVal = @(t,x) exp(-R*(p.Tmax-t)) * p.finalTimeVal(x);
        else
            p.boundaryVal = @(t,x) p.finalTimeVal(x);
        end
        switch p.bctype
        case 1
            p.DirichletBoundary = @(on_boundary, x) (x(on_boundary,1) < EPS | x(on_boundary,2) < EPS);
        case 2
            p.DirichletBoundary = @(on_boundary, x) (x(on_boundary,1) < EPS | x(on_boundary,2) < EPS | ( x(on_boundary,1) > p.xmax-EPS & x(on_boundary,2) > p.ymax-EPS));
        end

    case 'call'
        p.finalTimeVal = @(x) max(min(x,[],2)-F,0);
        p.boundaryVal  = @(t,x) zeros(n(x),1);
        p.DirichletBoundary = @(on_boundary, x) (x(on_boundary,1) < EPS | x(on_boundary,2));
    otherwise
        fprintf('Payoff type unknown')
    end

    switch p.bctype
    case 1
        p.diffusion{1, 1} = @(x) .5*sigmaV^2*x(:,1).^2;
        p.diffusion{1, 2} = @(x) .5*rho_VH *sigmaV*sigmaH*x(:,1).*x(:,2);
        p.diffusion{2, 1} = @(x) .5*rho_VH *sigmaV*sigmaH*x(:,1).*x(:,2);
        p.diffusion{2, 2} = @(x) .5*sigmaH^2*x(:,2).^2;
        p.convection{1} = @(x, t, u) R*x(:,1);
        p.convection{2} = @(x, t, u) R*x(:,2);
    case 2
        p.diffusion{1, 1} = @(x) .5*sigmaV^2*x(:,1).^2 .*(x(:,1) < p.xmax);
        p.diffusion{1, 2} = @(x) .5*rho_VH*sigmaV*sigmaH*x(:,1).*x(:,2) .* (x(:,1) < p.xmax) .* (x(:,2) < p.ymax);
        p.diffusion{2, 1} = @(x) .5*rho_VH*sigmaV*sigmaH*x(:,1).*x(:,2).* (x(:,1) < p.xmax) .* (x(:,2) < p.ymax);
        p.diffusion{2, 2} = @(x) .5*sigmaH^2*x(:,2).^2.*(x(:,2) < p.ymax);
        p.convection{1} = @(x, t, u) R*x(:,1) .* (x(:,1) < p.xmax);
        p.convection{2} = @(x, t, u) R*x(:,2) .* (x(:,2) < p.ymax);
    end

    switch p.type
    case 'European'
        solution = example_WorstAssetOption_Solution(F, R, sigmaH, sigmaV, rho_VH, p.Tmax, p.payoff);
        p.solution = @(x) solution(p.Tmin, x);
	    p.trivialControl = 1;
    case 'American'
        p.determineOptimalControl = @example_AmericanWorstAssetPut_Control;
	    p.trivialControl = 0;
    end

    p.control_dim = 1;
end
