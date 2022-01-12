function plotLabel(p, a, label)

if p.dim == 2
z = zeros(size(a.X));
z(label) = 1;
surf(a.x, a.y, reshape(z, p.nx+1, p.ny+1))
title('Boundary Label')
xlabel('x')
ylabel('y')
end % if p.dim == 2

if p.dim == 3
	error('to be implemented')
end % if p.dim == 3
