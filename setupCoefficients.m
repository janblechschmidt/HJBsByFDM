% This file provides the coefficients for the energy storage problem
function c = setupCoefficients(p)
EPS = 1e-10;
n = @(x) size(x,1);

if strcmp( p.mode, 'MinimumArrivalTime' )

	c.DirichletBC = 1;
	c.alpha = p.alpha;   % Coefficient for L2 cost term
	c.beta = p.beta;    % Coefficient for L1 cost term
	valx = 0.5*p.sigmax^2;
	valy = 0.5*p.sigmay^2;
	valz = 0.5*p.sigmaz^2;
	valxz = 0.5 * p.corrxz * p.sigmax * p.sigmaz;
	valxy = 0.5 * p.corrxy * p.sigmax * p.sigmay;
	valyz = 0.5 * p.corryz * p.sigmay * p.sigmaz;
	
	if p.dim == 2
		c.diffusion    = {@(x) valx*ones(n(x),1);
		                  @(x) valxy*ones(n(x),1);
		                  @(x) valxy*ones(n(x),1);
		                  @(x) valy*ones(n(x),1)};
	
		c.convection   = {@(x,t,u) +u + p.mux *ones(n(x),1);
		                  @(x,t,u) p.muy*ones(n(x),1)};
	else
		if p.dim == 3
			c.diffusion    = {@(x) valx*ones(n(x),1);
			                  @(x) valxy*ones(n(x),1);
			                  @(x) valxz*ones(n(x),1);
			                  @(x) valxy*ones(n(x),1);
			                  @(x) valy*ones(n(x),1);
			                  @(x) valyz*ones(n(x),1);
			                  @(x) valxz*ones(n(x),1);
			                  @(x) valyz*ones(n(x),1);
			                  @(x) valz*ones(n(x),1)};
	
		c.convection   = {@(x,t,u) +u + p.mux *ones(n(x),1);
		                  @(x,t,u) p.muy*ones(n(x),1);
		                  @(x,t,u) p.muz*ones(n(x),1)};
			% TODO: This is currently not used

		else
			error('Coefficients for dimension %i not implemented\n', p.dim)
		end % if p.dim == 3
	end % if p.dim == 2
	c.potential    = @(x) zeros(n(x),1);
	c.finalTimeVal = @(x) zeros(n(x),1);
	c.boundaryVal  = @(x) zeros(n(x),1);
	c.f            = @(x,u) (1 + 0.5 * c.alpha * u.^2 + c.beta * abs(u));

end % if strcmp( p.mode, 'MinimumArrivalTime' )

if strcmp( p.mode, 'EnergyStorage')

	c.DirichletBC = 0;

	if p.dim == 2

		if p.seasonality
			% c.meanPrice = @(t) 29.0 - 5.0*t/365 + 2.5*cos(2*pi*(1+t/365));
			c.meanPrice = @(t) 40.0 -3.0*t/365 + 2.5*cos(2*pi*(1+t/365));
		else 
			% c.meanPrice = @(t) 29.0;
			c.meanPrice = @(t) 40.0;
		end % if p.seasonality

		c.diffusion    = {@(x) 0.5*p.sigmap^2*ones(n(x),1);
											@(x) zeros(n(x),1);
											@(x) zeros(n(x),1);
											@(x) p.viscosity*ones(n(x),1)};
		c.convection   = {@(x,t,u) -p.kappap * (c.meanPrice(t) - x(:,1));
											@(x,t,u) -u};
		c.potential    = @(x) p.r * ones(n(x),1);
		c.finalTimeVal = @(x) zeros(n(x),1);
		c.boundaryVal  = @(x) zeros(n(x),1);
		c.alpha        = @(x) - p.K1 * sqrt(x(:,2));
		c.beta         = @(x) p.K2 * sqrt(1.0./(x(:,2) + p.K3) - 1.0 / p.K4);
		c.f2           = @(x) - p.coststorage  * x(:,2);
		c.f1           = @(x, u) (u > EPS) .* (- u .* x(:,1) * (1 + p.costpbuy) - p.costfbuy)  + (u < -EPS) .* (- u .* x(:,1) * (1 - p.costpsell) - p.costfsell);
		c.f            = @(x,u) (c.f1(x,u) + c.f2(x));
	
	else
		
		p.lambda = 0.5 * p.correlationcp * p.sigmap * p.sigmac;

		if p.seasonality

			c.meanPrice = @(t) 40.0 -3.0*t/365 + 2.5*cos(2*pi*(1+t/365));

			c.meanConsumption = @(t) 8000.0 - 100.0*t/365 + 1400.0*cos(5.93+2*pi*t/365) + 100*cos(5.52+4*pi*t/365);

		else 
			c.meanPrice = @(t) 40.0*ones(size(t));

			c.meanConsumption = @(t) 8000.0*ones(size(t));
		end % if p.seasonality

		c.diffusion    = {@(x) 0.5*p.sigmap^2*ones(n(x),1);
											@(x) zeros(n(x),1);
											@(x) p.lambda*ones(n(x),1);
											@(x) zeros(n(x),1);
											@(x) p.viscosity*ones(n(x),1);
											@(x) zeros(n(x),1);
											@(x) p.lambda*ones(n(x),1);
											@(x) zeros(n(x),1);
											@(x) 0.5*p.sigmac^2*ones(n(x),1)};
		c.convection   = {@(x,t,u) -p.kappap * (c.meanPrice(t) - x(:,1));
											@(x,t,u) -(u - x(:,3));
											@(x,t,u) -p.kappac * (c.meanConsumption(t) - x(:,3))};
		c.potential    = @(x) p.r * ones(n(x),1);
		% c.fillLevel = 0;
		c.fillLevel = 110000;
		c.penaltyFactor = 1.2;
		c.finalTimeVal = @(x) c.penaltyFactor * x(:,1).*(x(:,2)-c.fillLevel).*(x(:,2)-c.fillLevel < 0);
		c.boundaryVal  = @(x) zeros(n(x),1);
		c.alpha        = @(x) - p.K1 * sqrt(x(:,2)) + x(:,3);
		c.beta         = @(x) p.K2 * sqrt(1.0./(x(:,2) + p.K3) - 1.0 / p.K4) + x(:,3);
		c.f2           = @(x) - p.coststorage  * x(:,2);
		c.f1           = @(x, u) (u > EPS)   .* (- u      .* x(:,1) * (1 + p.costpbuy) - p.costfbuy)  + (u < -EPS) .* (- u .* x(:,1) * (1 - p.costpsell) - p.costfsell);
		c.f0           = @(x) (x(:,3) > EPS) .* (- x(:,3) .* x(:,1) * (1 + p.costpbuy) - p.costfbuy);
		if p.opportunity
			c.f            = @(x,u) c.f1(x,u) + c.f2(x) - c.f0(x);
		else
			c.f            = @(x,u) c.f1(x,u) + c.f2(x);
		end
	
	end % if p.dim == 2
end % if strcmp( p.mode, 'EnergyStorage')
end % function setupCoefficients(p)
