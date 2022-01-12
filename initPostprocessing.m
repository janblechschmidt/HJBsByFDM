function p = initPostprocessing(p, a)

switch p.dim
case 2
	% Obtain triangulation of quadrilateral qrid
	p.tri = delaunay(a.X, a.Y);
	p.Tri = delaunayTriangulation(a.X, a.Y);
	p.numtri = size(p.tri,1);
	% Ensure correct dimensions
	if p.numtri ~= p.nx*p.ny*2
		error('Wrong size of triangulation')
	end
case 3
	% Obtain triangulation of quadrilateral qrid
	p.tri = delaunay(a.X, a.Y, a.Z);
	p.Tri = delaunayTriangulation(a.X, a.Y, a.Z);
	p.numtri = size(p.tri,1);
	% Ensure correct dimensions
	if p.numtri ~= p.nx*p.ny*p.nz*6
		error('Wrong size of triangulation')
	end
end % switch p.dim
% triplot(tri, a.X, a.Y);
