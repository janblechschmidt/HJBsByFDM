function plotConvection(p, a, t, u)

switch p.dim
case 2
	myplot = @(fv) surf(a.x, a.y, reshape(fv, p.nx+1, p.ny+1),'FaceColor','interp');
case 3
	[X,Y,Z] = meshgrid(a.x, a.y, a.z);
	myplot = @(fv) slice(X, Y, Z, p.vec2Mat(fv),.5*(p.xmax-p.xmin),.5*(p.ymax-p.ymin),[.2 .8]*(p.zmax-p.zmin));
end

% Plot convection
for i = 1:length(p.convection)
	fh = p.convection{i};
	subplot(1, p.dim,i);
	fv = fh(a.XY, t, u);
	myplot(fv);
	mfv = mean(fv);
	colormap hsv
	caxis([min(min(fv),0),max(fv)]);

	ts = sprintf('Convection %d (Mean: %5.2f)', i, mfv);
	title(ts);
	xlabel('x');
	ylabel('y');
end

end % function plotCoefficients(p, c, a)
