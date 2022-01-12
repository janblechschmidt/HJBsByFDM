function p = example_EnergyStorage_Problem(p)

EPS = 1e-10;

if ~isfield(p, 'nx')
    p.nx = 20;
end
if ~isfield(p, 'ny')
    p.ny = 20;
end
if ~isfield(p, 'nz')
    p.nz = 20;
end
if ~isfield(p, 'nt')
    p.nt = 365;
end
if ~isfield(p, 'Tmin')
    p.Tmin = 0;
end
if ~isfield(p, 'Tmax')
    p.Tmax = 365;
end
p.pmin = 0;
p.pmax = 80;
p.qmin = 0;
p.qmax = 550000;
p.cmin = 0;
p.cmax = 16000;
p.qfactor = 1;
p.qshape = 0.25;
p.qout = -27000;
p.qin = 23000;
% p.costpsell = 0;
% p.costcsell = 0;
% p.costfsell = 0;
% p.costpbuy = 0;
% p.costcbuy = 0;
% p.costfbuy = 0;
% p.coststorage = 0;
p.costpsell = 0.05;
p.costcsell = 1e4;
p.costfsell = 0;
p.costpbuy = 0.05;
p.costcbuy = 1e4;
p.costfbuy = 0;
p.coststorage = 0.001;
if ~isfield(p, 'seasonality')
    p.seasonality = 1;
end
if ~isfield(p, 'opportunity')
    p.opportunity = 1;
end
p.sigmap = 0.5;
p.sigmac = 100;
if ~isfield(p, 'correlationcp')
    p.correlationcp = 0.3;
end
p.kappap = 0.02;
p.kappac = 0.007;
p.r = 0;

p.xmax = p.pmax;
p.xmin = p.pmin;
p.ymax = p.qmax;
p.ymin = p.qmin;

p.control_dim = 1;

if p.dim == 3
    p.zmax = p.cmax;
    p.zmin = p.cmin;
end

% Calculate some derived parameters for energy storage itself
p.K1 = - p.qout / sqrt(p.qmax);
p.K3 = p.qmax * p.qshape;
p.K4 = p.qmax + p.K3;
p.K2 = p.qin / sqrt(1.0/p.K3 - 1.0/p.K4);

% No Dirichlet boundary conditions
n = @(x) size(x,1);

p.DirichletBC = 0;

% Potential term
p.potential    = @(x) p.r * ones(n(x),1);

% Final conditions
p.fillLevel = 110000;
if ~isfield(p, 'bctype')
    p.bctype = 3;
end

switch p.bctype
    case 1
        p.finalTimeVal = @(x) zeros(n(x),1);
    case 2
        p.penaltyFactor = .9;
        p.finalTimeVal = @(x) p.penaltyFactor * x(:,1).*x(:,2);
    case 3 
        p.penaltyFactor = 2.;
        p.finalTimeVal = @(x) p.penaltyFactor * x(:,1).*(x(:,2)-p.fillLevel).*(x(:,2)-p.fillLevel < 0);
end
p.boundaryVal  = @(x) zeros(n(x),1);

% Running gain/cost
p.f1 = @(x, u) (u > EPS) .* (- u .* ( x(:,1) * (1 + p.costpbuy) + p.costfbuy) - p.costcbuy ) ...
    + (u < -EPS) .* (- u .* (x(:,1) * (1 - p.costpsell) - p.costfsell) - p.costcsell);
p.f2           = @(x) - p.coststorage  * x(:,2);

p.deterministicDimension = [2];

if ~isfield(p, 'meanPconst')
    p.meanPconst = 40.0;
end
if ~isfield(p, 'meanCconst')
    p.meanCconst = 8000.0;
end

if p.dim == 2

	if p.seasonality
		% p.meanPrice = @(t) 29.0 - 5.0*t/365 + 2.5*cos(2*pi*(1+t/365));
		p.meanPrice = @(t) p.meanPconst -3.0*t/365 + 2.5*cos(2*pi*(1+t/365));
	else 
		% p.meanPrice = @(t) 29.0;
		p.meanPrice = @(t) p.meanPconst;
	end % if p.seasonality

    % p.diffusion{1,1} = @(x) 0.5*p.sigmap^2*ones(n(x),1) .* (x(:,1) < p.xmax-EPS) .* (x(:,1)>p.xmin+EPS);
    % p.diffusion{1,1} = @(x) 0.5*p.sigmap^2*ones(n(x),1) .* (x(:,1) < p.xmax-EPS);
    p.diffusion{1,1} = @(x) 0.5*p.sigmap^2*ones(n(x),1);
    p.diffusion{1,2} = @(x) zeros(n(x),1);
    p.diffusion{2,1} = @(x) zeros(n(x),1);
    p.diffusion{2,2} = @(x) zeros(n(x),1);

    p.convection{1} = @(x,t,u) p.kappap * (p.meanPrice(t) - x(:,1));
    p.convection{2} = @(x,t,u) u;

	p.alpha        = @(x) - p.K1 * sqrt(x(:,2));
	p.beta         = @(x) p.K2 * sqrt(1.0./(x(:,2) + p.K3) - 1.0 / p.K4);
	p.f            = @(x,u) (p.f1(x,u) + p.f2(x));

    p.qoi = [40, 110000];

else
	
	p.lambda = 0.5 * p.correlationcp * p.sigmap * p.sigmac;

	if p.seasonality
		p.meanPrice = @(t) p.meanPconst -3.0*t/365 + 2.5*cos(2*pi*(1+t/365));
		p.meanConsumption = @(t) p.meanCconst - 100.0*t/365 + 1400.0*cos(5.93+2*pi*t/365) + 100*cos(5.52+4*pi*t/365);
        p.timeConstantConvection = 0;

	else 
		p.meanPrice = @(t) p.meanPconst * ones(size(t));
		p.meanConsumption = @(t) p.meanCconst * ones(size(t));
        p.timeConstantConvection = 1;

	end % if p.seasonality

    p.diffusion{1,1} = @(x) 0.5*p.sigmap^2*ones(n(x),1);
    p.diffusion{1,2} = @(x) zeros(n(x),1);
    p.diffusion{1,3} = @(x) p.lambda*ones(n(x),1);
    p.diffusion{2,1} = @(x) zeros(n(x),1);
    p.diffusion{2,2} = @(x) zeros(n(x),1);
    p.diffusion{2,3} = @(x) zeros(n(x),1);
    p.diffusion{3,1} = @(x) p.lambda*ones(n(x),1);
    p.diffusion{3,2} = @(x) zeros(n(x),1);
    p.diffusion{3,3} = @(x) 0.5*p.sigmac^2*ones(n(x),1);

    p.convection{1} = @(x,t,u) p.kappap * (p.meanPrice(t) - x(:,1));
    p.convection{2} = @(x,t,u) u - x(:,3);
    p.convection{3} = @(x,t,u) p.kappac * (p.meanConsumption(t) - x(:,3));

	p.alpha        = @(x) - p.K1 * sqrt(x(:,2)) + x(:,3);
	p.beta         = @(x) p.K2 * sqrt(1.0./(x(:,2) + p.K3) - 1.0 / p.K4) + x(:,3);

	% p.f0           = @(x) (x(:,3) > EPS) .* (- x(:,3) .* x(:,1) * (1 + p.costpbuy) - p.costfbuy);
	p.f0           = @(x) p.f1(x, x(:,3));
	if p.opportunity
		p.f            = @(x,u) p.f1(x, u) + p.f2(x) - p.f0(x);
	else
		p.f            = @(x,u) p.f1(x, u) + p.f2(x);
    end

    p.qoi = [40, 110000, 8000];

    
    p.cntrls = @(x) [p.alpha(x), p.beta(x), zeros(size(x,1),1) + p.alpha(x) .* (p.alpha(x)>0)];

end % if p.dim == 2

p.determineOptimalControl = @example_EnergyStorage_Control;
