function [a,p] = setupStructure(p)

grow = @(x) 0.5*([x,0] + [0,x]);

refinement_level = p.refinement.level;

p.ht = (p.Tmax - p.Tmin) / p.nt;

if p.dim == 1

	% Setup spatial grid
	a.x = linspace(p.xmin, p.xmax, p.nx + 1);
    
	if refinement_level > 0
		fprintf('Refine cells\n')
        if ~isfield(p.refinement,'xlow')
            xlow = floor(1/3*(p.nx+1));
            xup  = ceil(2/3*(p.nx+1));
    		[a.x, p.nx] = refineCells(a.x, xlow:xup);
        else
            for i = 1:length(p.refinement.xlow)
                xlowi = p.refinement.xlow(i);
                xupi = p.refinement.xup(i);
                if p.refinement.mode == 'abs'
                    xlow = min(find(a.x>=xlowi));
                    xup = max(find(a.x<=xupi));
                else
                    xlow = floor(xlowi*(p.nx+1));
                    xup = ceil(xupi*(p.nx+1));
                end
    		    [a.x, p.nx] = refineCells(a.x, xlow:xup);
            end
        end
	end % if refinement_level > 0

    X = a.x;
    a.X = X(:);
    a.XY = [a.X];

	% Setup index vectors
	ix = 1:1:p.nx + 1;
	a.ix = ix(:);

	% Setup inner indices
	a.ixinner = find(a.ix ~= 1 & a.ix ~= p.nx + 1);
	
	% Setup boundaries
	a.ixlow = find(a.ix == 1);
	a.ixup = find(a.ix == p.nx + 1);

	a.iouter = find(a.ix == 1 | a.ix == p.nx + 1);
	a.iinner = find(a.ix ~= 1 & a.ix ~= p.nx + 1);

	% TODO: p.nx etc. should be replaced by a.nx
	a.nx = p.nx;

	% Set up spatial differences for differential operator
	a.dx = diff(a.x);
    tmpx = [a.dx, nan];

	a.Dx = tmpx(:);

	% Set up 'volume of cell' which is necessary for correct set up of L2 norm
	tmpx = grow(a.dx);
	a.H = tmpx;
	a.Hv = a.H(:);

end

if p.dim == 2

	% Setup spatial grid
	a.x = linspace(p.xmin, p.xmax, p.nx + 1);
	a.y = linspace(p.ymin, p.ymax, p.ny + 1);


	if refinement_level > 0
		fprintf('Refine cells\n')
        if ~isfield(p.refinement,'xlow')
            xlow = floor(1/3*(p.nx+1));
            xup  = ceil(2/3*(p.nx+1));
    		[a.x, p.nx] = refineCells(a.x, xlow:xup);
        else
            for i = 1:length(p.refinement.xlow)
                xlowi = p.refinement.xlow(i);
                xupi = p.refinement.xup(i);
                if p.refinement.mode == 'abs'
                    xlow = min(find(a.x>=xlowi));
                    xup = max(find(a.x<=xupi));
                else
                    xlow = floor(xlowi*(p.nx+1));
                    xup = ceil(xupi*(p.nx+1));
                end
    		    [a.x, p.nx] = refineCells(a.x, xlow:xup);
            end
        end

        if ~isfield(p.refinement,'ylow')
            ylow = floor(1/3*(p.ny+1));
            yup  = ceil(2/3*(p.ny+1));
    		[a.y, p.ny] = refineCells(a.y, ylow:yup);
        else
            for i = 1:length(p.refinement.ylow)
                ylowi = p.refinement.ylow(i);
                yupi = p.refinement.yup(i);
                if p.refinement.mode == 'abs'
                    ylow = min(find(a.y>=ylowi));
                    yup = max(find(a.y<=yupi));
                else
                    ylow = floor(ylowi*(p.ny+1));
                    yup = ceil(yupi*(p.ny+1));
                end
    		    [a.y, p.ny] = refineCells(a.y, ylow:yup);
            end
        end
	end % if refinement_level > 0

	% f = 1.5;
	% meshtrans = @(x) min(x) + (x/(max(x)-min(x))).^f * (max(x)-min(x));
	% a.x = meshtrans(a.x);
	% a.y = meshtrans(a.y);

	[X, Y] = meshgrid(a.x, a.y);
	a.X = X(:);
	a.Y = Y(:);
	a.XY = [a.X, a.Y];
	
	% Setup index vectors
	[ix, iy] = meshgrid([1:1:p.nx + 1], [1:1:p.ny + 1]);
	a.ix = ix(:);
	a.iy = iy(:);

	% Setup inner indices
	a.ixinner = find(a.ix ~= 1 & a.ix ~= p.nx + 1);
	a.iyinner = find(a.iy ~= 1 & a.iy ~= p.ny + 1);
	
	% Setup boundaries
	a.ixlow = find(a.ix == 1);
	a.ixup = find(a.ix == p.nx + 1);
	a.iylow = find(a.iy == 1);
	a.iyup = find(a.iy == p.ny + 1);

	% Setup inner of boundaries
	a.ixlowinner = find(a.ix == 1 & a.iy ~= 1 & a.iy ~= p.ny + 1) ;
	a.ixupinner  = find(a.ix == p.nx + 1 & a.iy ~= 1 & a.iy ~= p.ny + 1) ;
	a.iylowinner = find(a.iy == 1 & a.ix ~= 1 & a.ix ~= p.nx + 1);
	a.iyupinner  = find(a.iy == p.ny + 1 & a.ix ~= 1 & a.ix ~= p.nx + 1);

	a.iouter = find(a.ix == 1 | a.ix == p.nx + 1 | a.iy == 1 | a.iy == p.ny + 1);
	a.iinner = find(a.ix ~= 1 & a.ix ~= p.nx + 1 & a.iy ~= 1 & a.iy ~= p.ny + 1);

	% TODO: p.nx etc. should be replaced by a.nx
	a.nx = p.nx;
	a.ny = p.ny;

	% Set up spatial differences for differential operator
	a.dx = diff(a.x);
	a.dy = diff(a.y);
	[tmpx, tmpy] = meshgrid([a.dx,nan], [a.dy,nan]);

	a.Dx = tmpx(:);
	a.Dy = tmpy(:);

	% Set up 'volume of cell' which is necessary for correct set up of L2 norm
	[tmpx, tmpy] = meshgrid(grow(a.dx), grow(a.dy));
	a.H = tmpx .* tmpy;
	a.Hv = a.H(:);

end % if p.dim == 2

if p.dim == 3
	% Setup spatial grid
	a.x = linspace(p.xmin, p.xmax, p.nx + 1);
	a.y = linspace(p.ymin, p.ymax, p.ny + 1);
	a.z = linspace(p.zmin, p.zmax, p.nz + 1);

	if refinement_level > 0

		fprintf('Refine cells\n')

        if ~isfield(p.refinement,'xlow')
            xlow = floor(1/3*(p.nx+1));
            xup  = ceil(2/3*(p.nx+1));
    		[a.x, p.nx] = refineCells(a.x, xlow:xup);
        else
            for i = 1:length(p.refinement.xlow)
                xlowi = p.refinement.xlow(i);
                xupi = p.refinement.xup(i);
                if p.refinement.mode == 'abs'
                    xlow = min(find(a.x>=xlowi));
                    xup = max(find(a.x<=xupi));
                else
                    xlow = floor(xlowi*(p.nx+1));
                    xup = ceil(xupi*(p.nx+1));
                end
    		    [a.x, p.nx] = refineCells(a.x, xlow:xup);
            end
        end

        if ~isfield(p.refinement,'ylow')
            ylow = floor(1/3*(p.ny+1));
            yup  = ceil(2/3*(p.ny+1));
    		[a.y, p.ny] = refineCells(a.y, ylow:yup);
        else
            for i = 1:length(p.refinement.ylow)
                ylowi = p.refinement.ylow(i);
                yupi = p.refinement.yup(i);
                if p.refinement.mode == 'abs'
                    ylow = min(find(a.y>=ylowi));
                    yup = max(find(a.y<=yupi));
                else
                    ylow = floor(ylowi*(p.ny+1));
                    yup = ceil(yupi*(p.ny+1));
                end
    		    [a.y, p.ny] = refineCells(a.y, ylow:yup);
            end
        end

        if ~isfield(p.refinement,'zlow')
            zlow = floor(1/3*(p.nz+1));
            zup  = ceil(2/3*(p.nz+1));
    		[a.z, p.nz] = refineCells(a.z, zlow:zup);
        else
            for i = 1:length(p.refinement.zlow)
                zlowi = p.refinement.zlow(i);
                zupi = p.refinement.zup(i);
                if p.refinement.mode == 'abs'
                    zlow = min(find(a.z>=zlowi));
                    zup = max(find(a.z<=zupi));
                else
                    zlow = floor(zlowi*(p.nz+1));
                    zup = ceil(zupi*(p.nz+1));
                end
    		    [a.z, p.nz] = refineCells(a.z, zlow:zup);
            end
        end

	end % if refinement_level > 0

	% f = 1.5;
	% meshtrans = @(x) min(x) + (x/(max(x)-min(x))).^f * (max(x)-min(x));
	% a.x = meshtrans(a.x);
	% a.y = meshtrans(a.y);
	% a.z = meshtrans(a.z);

	[X, Y, Z] = meshgrid(a.x, a.y, a.z);
	a.X = X(:);
	a.Y = Y(:);
	a.Z = Z(:);
	a.XY = [a.X, a.Y, a.Z];

	% Setup index vectors
	[ix, iy, iz] = meshgrid([1:1:p.nx + 1], [1:1:p.ny + 1], [1:1:p.nz + 1]);
	a.ix = ix(:);
	a.iy = iy(:);
	a.iz = iz(:);

	% Setup inner indices
	a.ixinner = find(a.ix ~= 1 & a.ix ~= p.nx + 1);
	a.iyinner = find(a.iy ~= 1 & a.iy ~= p.ny + 1);
	a.izinner = find(a.iz ~= 1 & a.iz ~= p.nz + 1);

	% TODO: A lot of the following indices is not really needed.
	
	% Setup faces
	a.ixlow = find(a.ix == 1);
	a.ixup  = find(a.ix == p.nx + 1);
	a.iylow = find(a.iy == 1);
	a.iyup  = find(a.iy == p.ny + 1);
	a.izlow = find(a.iz == 1);
	a.izup  = find(a.iz == p.nz + 1);

	% Setup inner of faces
	a.ixlowinner = find(a.ix == 1 & a.iy ~= 1 & a.iy ~= p.ny + 1 & a.iz ~= 1 & a.iz ~= p.nz + 1);
	a.ixupinner  = find(a.ix == p.nx + 1 & a.iy ~= 1 & a.iy ~= p.ny + 1 & a.iz ~= 1 & a.iz ~= p.nz + 1);
	a.iylowinner = find(a.iy == 1 & a.ix ~= 1 & a.ix ~= p.nx + 1 & a.iz ~= 1 & a.iz ~= p.nz + 1);
	a.iyupinner  = find(a.iy == p.ny + 1 & a.ix ~= 1 & a.ix ~= p.nx + 1 & a.iz ~= 1 & a.iz ~= p.nz + 1);
	a.izlowinner = find(a.iz == 1 & a.ix ~= 1 & a.ix ~= p.nx + 1 & a.iy ~= 1 & a.iy ~= p.ny + 1);
	a.izupinner  = find(a.iz == p.nz + 1 & a.ix ~= 1 & a.ix ~= p.nx + 1 & a.iy ~= 1 & a.iy ~= p.ny + 1);

	% Setup edges
	a.ixlow_iylow = find(a.ix == 1 & a.iy == 1);
	a.ixlow_iyup = find(a.ix == 1 & a.iy == p.ny + 1);
	a.ixup_iylow  = find(a.ix == p.nx + 1 & a.iy == 1);
	a.ixup_iyup  = find(a.ix == p.nx + 1 & a.iy == p.ny + 1);
	
	a.ixlow_izlow = find(a.ix == 1 & a.iz == 1);
	a.ixlow_izup = find(a.ix == 1 & a.iz == p.nz + 1);
	a.ixup_izlow  = find(a.ix == p.nx + 1 & a.iz == 1);
	a.ixup_izup  = find(a.ix == p.nx + 1 & a.iz == p.nz + 1);

	a.iylow_izlow = find(a.iy == 1 & a.iz == 1);
	a.iylow_izup = find(a.iy == 1 & a.iz == p.nz + 1);
	a.iyup_izlow  = find(a.iy == p.ny + 1 & a.iz == 1);
	a.iyup_izup  = find(a.iy == p.ny + 1 & a.iz == p.nz + 1);

	
	% Setup inner of edges
	a.ixlow_iylow_inner = find(a.ix == 1 & a.iy == 1 & a.iz ~= 1 & a.iz ~= p.nz + 1);
	a.ixlow_iyup_inner = find(a.ix == 1 & a.iy == p.ny + 1 & a.iz ~= 1 & a.iz ~= p.nz + 1);
	a.ixup_iylow_inner  = find(a.ix == p.nx + 1 & a.iy == 1 & a.iz ~= 1 & a.iz ~= p.nz + 1);
	a.ixup_iyup_inner  = find(a.ix == p.nx + 1 & a.iy == p.ny + 1 & a.iz ~= 1 & a.iz ~= p.nz + 1);
	
	a.ixlow_izlow_inner = find(a.ix == 1 & a.iz == 1 & a.iy ~= 1 & a.iy ~= p.ny + 1);
	a.ixlow_izup_inner = find(a.ix == 1 & a.iz == p.nz + 1 & a.iy ~= 1 & a.iy ~= p.ny + 1);
	a.ixup_izlow_inner  = find(a.ix == p.nx + 1 & a.iz == 1 & a.iy ~= 1 & a.iy ~= p.ny + 1);
	a.ixup_izup_inner  = find(a.ix == p.nx + 1 & a.iz == p.nz + 1 & a.iy ~= 1 & a.iy ~= p.ny + 1);

	a.iylow_izlow_inner = find(a.iy == 1 & a.iz == 1 & a.ix ~= 1 & a.ix ~= p.nx + 1);
	a.iylow_izup_inner = find(a.iy == 1 & a.iz == p.nz + 1 & a.ix ~= 1 & a.ix ~= p.nx + 1);
	a.iyup_izlow_inner  = find(a.iy == p.ny + 1 & a.iz == 1 & a.ix ~= 1 & a.ix ~= p.nx + 1);
	a.iyup_izup_inner  = find(a.iy == p.ny + 1 & a.iz == p.nz + 1 & a.ix ~= 1 & a.ix ~= p.nx + 1);
	
	a.iouter = find(a.ix == 1 | a.ix == p.nx + 1 | a.iy == 1 | a.iy == p.ny + 1 | a.iz == 1 | a.iz == p.nz + 1);
	a.iinner = find(a.ix ~= 1 & a.ix ~= p.nx + 1 & a.iy ~= 1 & a.iy ~= p.ny + 1 & a.iz ~= 1 & a.iz ~= p.nz + 1);

	
	a.dx = diff(a.x);
	a.dy = diff(a.y);
	a.dz = diff(a.z);

	[tmpx, tmpy, tmpz] = meshgrid([a.dx,nan], [a.dy,nan], [a.dz,nan]);
	a.Dx = tmpx(:);
	a.Dy = tmpy(:);
	a.Dz = tmpz(:);

	a.nx = p.nx;
	a.ny = p.ny;
	a.nz = p.nz;

	[tmpx, tmpy, tmpz] = meshgrid(grow(a.dx), grow(a.dy), grow(a.dz));

	a.H = tmpx .* tmpy .* tmpz;
	a.Hv = a.H(:);

end % if p.dim == 3

switch p.dim
case 1
    p.ndofs = p.nx + 1;
	p.vec2Mat = @(v) v;
    a.get_qoi_val = @(w) interp1(a.X, w, p.qoi);

case 2
	p.ndofs = (p.nx+1)*(p.ny+1);
	p.vec2Mat = @(v) reshape(v,[p.ny+1,p.nx+1]);
    a.get_qoi_val = @(w) interp2(p.vec2Mat(a.X),p.vec2Mat(a.Y),p.vec2Mat(w),p.qoi(1),p.qoi(2) );

case 3
	p.ndofs = (p.nx+1)*(p.ny+1)*(p.nz+1);
	p.vec2Mat = @(v) reshape(v,[p.ny+1,p.nx+1,p.nz+1]);
    a.get_qoi_val = @(w) interp3(p.vec2Mat(a.X), p.vec2Mat(a.Y), p.vec2Mat(a.Z), p.vec2Mat(w), p.qoi(1),p.qoi(2), p.qoi(3));
end % switch p.dim

% Dirichlet boundary indices
if p.DirichletBC == 1
    a.bc_idxs = p.DirichletBoundary(a.iouter, a.XY);
    a.bc_idxs = a.iouter(find(a.bc_idxs));
    a.dof_idxs = setdiff(1:p.ndofs,a.bc_idxs);
else
    a.bc_idxs = [];
    a.dof_idxs = 1:p.ndofs;
end


% Setup time grid
a.t = linspace(p.Tmin, p.Tmax, p.nt + 1);
a.iall   = (1:p.ndofs)'; 
