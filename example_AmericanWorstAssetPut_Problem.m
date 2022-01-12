function p = example_AmericanWorstAssetPut_Problem(p)
    
    n = @(x) size(x,1);
    EPS = 1e-12;

    p.nt = 320;
    p.Tmin = 0.5;
    p.Tmax = 1;
    p.nx = 100;
    p.xmin = 0;
    p.xmax = 100;
    p.ny = 100;
    p.ymin = 0;
    p.ymax = 100;

    % tau = 0.5;
    H = 40;
    V = 40;
    F = 40;
    R = 0.05;
    sigmaH = 0.3;
    sigmaV = 0.3;
    rho_VH = 0.5;
    p.qoi = [V, H];

    p.potential    = @(x, u) R * ones(n(x),1);
    p.f            = @(x, u) zeros(n(x),1);

    p.DirichletBC = 1;
    p.DirichletBoundary = @(on_boundary, x) (x(on_boundary,1) < EPS | x(on_boundary,2) < EPS);

    % Put
    p.boundaryVal = @(t,x) exp(-R*(p.Tmax-t))*max(F-min(x,[],2),0);
    p.finalTimeVal = @(x) max(F-min(x,[],2),0);
    
    % Call
    % p.finalTimeVal = @(x) max(min(x,[],2)-F,0);
    % p.boundaryVal  = @(x) zeros(n(x),1);


    p.diffusion{1, 1} = @(x,u) .5*sigmaV^2*x(:,1).^2;
    p.diffusion{1, 2} = @(x,u) .5*rho_VH*sigmaV*sigmaH*x(:,1).*x(:,2);
    p.diffusion{2, 1} = @(x,u) .5*rho_VH*sigmaV*sigmaH*x(:,1).*x(:,2);
    p.diffusion{2, 2} = @(x,u) .5*sigmaH^2*x(:,2).^2;

    p.convection{1} = @(x, t, u) R*x(:,1);
    p.convection{2} = @(x, t, u) R*x(:,2);

    p.control_dim = 1;

    p.determineOptimalControl = @example_AmericanWorstAssetPut_Control;
end
