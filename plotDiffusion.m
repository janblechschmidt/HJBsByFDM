function plotDiffusion(p, a)

switch p.dim
case 2
	myplot = @(fv) surf(a.x, a.y, reshape(fv, p.nx+1, p.ny+1));
case 3
	[X,Y,Z] = meshgrid(a.x, a.y, a.z);
	myplot = @(fv) slice(X, Y, Z, p.vec2Mat(fv),.5*(p.xmax-p.xmin),.5*(p.ymax-p.ymin),[.2 .8]*(p.zmax-p.zmin));
end

idx = 0;
for i = 1:p.dim
    for j = 1:p.dim
    	fh = p.diffusion{i, j};
        idx = idx + 1;
    	subplot(p.dim, p.dim, idx);
    	fv = fh(a.XY);
    	myplot(fv);
    	mfv = mean(fv);
    	colormap hsv
    	caxis([min(min(fv),0),max(fv)]);
    
    	ts = sprintf('Diffusion %d (Mean: %5.2f)', i, mfv);
    	title(ts);
    	xlabel = 'x';
    	ylable = 'y';
    end
end
