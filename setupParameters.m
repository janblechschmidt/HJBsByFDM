% This file sets all parameters for the energy storage problem
function p = setupParameters(p)

if strcmp( p.mode, 'MinimumArrivalTime' )
	% p = parseFEniCSxml('HJBMinimumArrivalTime.xml',p);
    p = setupProblem_MinimumArrivalTime(2)
else
	if strcmp( p.mode, 'EnergyStorage' )

		if p.dim == 2
			p = parseFEniCSxml('HJBEnergyStorage2d.xml',p);
		
		end % if p.dim == 2

		if p.dim == 3

			p = parseFEniCSxml('HJBEnergyStorage3d.xml',p);

		end % if p.dim == 3

		p.xmax = p.pmax;
		p.xmin = p.pmin;
		p.ymax = p.qmax;
		p.ymin = p.qmin;
	
		% Calculate some derived parameter for energy storage itself
		p.K1 = - p.qout / sqrt(p.qmax);
		p.K3 = p.qmax * p.qshape;
		p.K4 = p.qmax + p.K3;
		p.K2 = p.qin / sqrt(1.0/p.K3 - 1.0/p.K4);

		if p.dim == 3
			p.zmax = p.cmax;
			p.zmin = p.cmin;
		end

	end % if strcmp( p.mode, 'EnergyStorage' )
end % if strcmp( p.mode, 'MinimumArrivalTime' )
end % function p = setupParameters(p)
