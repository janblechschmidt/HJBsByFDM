function p = example_MinimumArrivalTime_Problem(p)

    if ~isfield(p, 'alpha')
        p.alpha = 0.5;
    end

    if ~isfield(p, 'beta')
        p.beta = 0.5;
    end

    n = @(x) size(x,1);

    p.potential    = @(x) zeros(n(x),1);
    p.finalTimeVal = @(x) zeros(n(x),1);
    p.DirichletBC = 1;
    p.DirichletBoundary = @(on_boundary, x) ones(size(on_boundary));
    p.boundaryVal  = @(t,x) zeros(n(x),1);
    p.f            = @(x,u) (1 + 0.5 * p.alpha * sum(u.^2, 2) + p.beta * sum(abs(u), 2));
    p.control_dim = p.dim;
    % Umin and Umax are only relevant if alpha = 0
    p.Umin = -1;
    p.Umax = 1;
    p.qoi = zeros(1,p.dim);

    if p.dim >= 1
        if ~isfield(p, 'nt')
            p.nt = 32;
        end
        p.Tmin = 0;
        p.Tmax = 1;

        if ~isfield(p, 'nx')
            p.nx = 1000;
        end
        if ~isfield(p, 'xmin')
            p.xmin = 0;
        end
        if ~isfield(p, 'xmax')
            p.xmax = 1;
        end

        if ~isfield(p, 'sigmax')
            p.sigmax = 0.5;
        end
        p.mux = 0.0;

        valxx = 0.5*p.sigmax^2;

        p.diffusion{1, 1} = @(x) valxx*ones(n(x),1);
        p.convection{1} = @(x, t, u) u(:, 1) + p.mux *ones(n(x),1);


    end % if p.dim >= 1

    if p.dim >= 2

        if ~isfield(p, 'nx')
            p.nx = 20;
        end
        if ~isfield(p, 'ny')
            p.ny = 20;
        end
        if ~isfield(p, 'ymin')
            p.ymin = 0;
        end
        if ~isfield(p, 'ymax')
            p.ymax = 1;
        end

        if ~isfield(p, 'sigmay')
            p.sigmay = 0.5;
        end

        if ~isfield(p, 'corrxy')
            p.corrxy = 0.0;
        end

        valyy = 0.5*p.sigmay^2;
        valxy = 0.5 * p.corrxy * p.sigmax * p.sigmay;

        p.diffusion{1, 2} = @(x) valxy*ones(n(x),1);
        p.diffusion{2, 1} = @(x) valxy*ones(n(x),1);
        p.diffusion{2, 2} = @(x) valyy*ones(n(x),1);

        % p.diffusion{1, 2} = @(x) (2*(sqrt(x(:,1).^2+x(:,2).^2)<.5)-1).*valxy.*ones(n(x),1);
        % p.diffusion{2, 1} = @(x) (2*(sqrt(x(:,1).^2+x(:,2).^2)<.5)-1).*valxy.*ones(n(x),1);

        p.muy = 0.0;
        p.convection{2} = @(x,t,u) u(:, 2) + p.muy*ones(n(x),1);
        
    end % if p.dim >= 2

    if p.dim >= 3

        if ~isfield(p, 'nx')
            p.nx = 10;
        end
        if ~isfield(p, 'ny')
            p.ny = 10;
        end
        if ~isfield(p, 'nz')
            p.nz = 10;
        end
        if ~isfield(p, 'zmin')
            p.zmin = 0;
        end
        if ~isfield(p, 'zmax')
            p.zmax = 1;
        end

        if ~isfield(p, 'sigmaz')
            p.sigmaz = 0.5;
        end
        p.muz = 0.0;

        p.corrxz = 0.0;
        p.corryz = 0.0;

        valzz = 0.5*p.sigmaz^2;
        valxz = 0.5 * p.corrxz * p.sigmax * p.sigmaz;
        valyz = 0.5 * p.corryz * p.sigmay * p.sigmaz;

        p.diffusion{1, 3} = @(x) valxz*ones(n(x),1);
        p.diffusion{2, 3} = @(x) valyz*ones(n(x),1);
        p.diffusion{3, 1} = @(x) valxz*ones(n(x),1);
        p.diffusion{3, 2} = @(x) valyz*ones(n(x),1);
        p.diffusion{3, 3} = @(x) valzz*ones(n(x),1);

        p.convection{3} = @(x,t,u) u(:, 3) + p.muz*ones(n(x),1);
        
    end % if p.dim >= 3

    p.determineOptimalControl = @example_MinimumArrivalTime_Control;
