function exportVTKSolution(p, a, v, u, name)

if nargin < 5
	name = 'solution.vtk';
end % if nargin < 5

if p.dim == 2
	% [X,Y,Z] = meshgrid(a.x, a.y, [0]);
	[X,Y,Z] = meshgrid(a.x/p.xmax, a.y/p.ymax, [0]);
	vtkwrite([p.outputdir '/' name], 'structured_grid', X, Y, Z,...
	'scalars', 'Value', p.vec2Mat(v), ...
	'scalars', 'Control', p.vec2Mat(u),'binary');
end % if p.dim == 2
if p.dim == 3
	% [X,Y,Z] = meshgrid(a.x, a.y, a.z);
	[X,Y,Z] = meshgrid(a.x/p.xmax, a.y/p.ymax, a.z/p.zmax);
	vtkwrite([p.outputdir '/' name], 'structured_grid', X, Y, Z, ...
	'scalars', 'Value', p.vec2Mat(v), ...
	'scalars', 'Control', p.vec2Mat(u),'binary');
end % if p.dim == 3
