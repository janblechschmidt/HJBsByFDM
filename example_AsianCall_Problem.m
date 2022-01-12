function p = example_AsianCall_Problem(p)

    n = @(x) size(x,1);
    EPS = 1e-12;

    % Set bounds of domain
    if ~isfield(p, 'nx')
        p.nx = 20;
    end
    p.xmin = 0;
    if ~isfield(p, 'xmax')
        p.xmax = 200;
    end
    if ~isfield(p, 'ny')
        p.ny = 20;
    end
    p.ymin = 0;
    if ~isfield(p, 'ymax')
        p.ymax = 200;
    end

    % Time horizon
    if ~isfield(p, 'nt')
        p.nt = 100;
    end
    % p.Tmax = 1.;
    % p.Tmin = 0.75;
    p.Tmax = 0.25;
    p.Tmin = 0.0;

    % Dirichlet boundary
    p.DirichletBC = 1;
    p.DirichletBoundary = @(onBoundary,x)  x(onBoundary,1) >= p.xmax-1e-5 & x(onBoundary,2) >= p.ymax - 1e-5;

    % Values from Zvan
    sigmaS = 0.1;
    r = 0.1;
    K = 100;

    % Zeroth order term
    p.potential    = @(x) r * ones(n(x),1);
    % .* (x(:,1)<p.xmax);

    p.deterministicDimension = 2;

    % First order term
    p.convection{1} = @(x, t, u) r*x(:,1) .* (x(:,1) > p.xmin) .* (x(:,1)<p.xmax);
    p.convection{2} = @(x, t, u) posConvection(t,x,(p.Tmax-p.Tmin)/p.nt);
    % @(t,x) (x(:,1) - x(:,2)) / (p.Tmax-t)};
    %@(t,x) (x(:,1) - x(:,2)) / t .* (t > 0)};

    % Second order term
    p.diffusion{1,1} = @(x) 0.5*sigmaS^2*x(:,1).^2 .* (x(:,1) > p.xmin)  .*(x(:,1)<p.xmax);
    p.diffusion{1,2} = @(x) zeros(n(x),1);
    p.diffusion{2,1} = @(x) zeros(n(x),1);
    p.diffusion{2,2} = @(x) zeros(n(x),1);

    p.boundaryVal = @(t,x) max(x(:,2)-exp(-r*(p.Tmax-t))*K,0);
    p.finalTimeVal = @(x) p.boundaryVal(p.Tmax,x);
    p.f = @(t,x) zeros(n(x),1);

    p.control_dim = 1;
    p.qoi = 100*ones(1,p.dim);
    p.timeConstantConvection = 0;
end

function ret = posConvection(t,x,dt)
    if t == 0
        t = dt;
        % ret = ones(size(x,1),1) .* ((x(:,1) > x(:,2)) - (x(:,1) < x(:,2)));
    end
    ret = (x(:,1) - x(:,2)) / t;
end
