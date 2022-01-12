function plotValue(p, a, w, lspace, suffix)
if nargin < 5
    suffix = '';
end
P = 20;
if p.dim == 1
	plot(a.x, p.vec2Mat(w));
	title('Value')
	xlabel('x')
	ylabel('y')

end % if p.dim == 1

if p.dim == 2

    % if p.writeValueFunction
    %     f = figure('visible','off');
    % end
	% azsurf = surf(a.x, a.y, p.vec2Mat(w),'EdgeAlpha',0.1,'LineWidth',0.001);
    % alpha(azsurf,0.9);
	% azsurf = surf(a.x, a.y, p.vec2Mat(w),'EdgeColor','none','LineWidth',0.001);
    % alpha(azsurf,0.9);
	azsurf = contour(a.x, a.y, p.vec2Mat(w),10,'LineWidth',1);
    caxis([0,1])
    colormap viridis;
    cb = colorbar('eastoutside');
    cb.Label.Interpreter='latex';
	% title('Value')
    % view([140,40])
    % axis([p.xmin, p.xmax, p.ymin, p.ymax,0,0.1])
    if isfield(p, 'xlabel')
        xlab = p.xlabel;
    else
        xlab = '$x^1$';
    end
	xlabel(xlab,'Interpreter','latex')
    if isfield(p, 'ylabel')
        ylab = p.ylabel;
    else
        ylab = '$x^2$';
    end
	ylabel(ylab,'Interpreter','latex')
    if p.writeValueFunction
        fname = strcat(p.outputdir,'/', p.prefix, sprintf('_value_theta_%02d_dofs_%05d',p.theta*10,length(a.dof_idxs)));
        if ~strcmp(suffix, '')
            fname = strcat(fname,'_',suffix);
        end
        set(gcf, 'Units', 'inches');
        set(gcf, 'PaperPosition', [0 0 P P]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
        print(gcf,fname,'-dpng','-r300'); 
    end
    % saveas(gcf, 'test.pdf') %Save figure
    % system('pdfcrop test.pdf test.pdf');

    % Matlab2tikz is too slow to compile
    % and needs too much memory
    % matlab2tikz('tmp.tex')

end % if p.dim == 2

if p.dim == 3

	if nargin < 5
		% lspace = round(p.nz/2);
		lspace = [round(p.nz/4), round(p.nz/2), round(p.nz*3/4)];
    end

    if p.writeValueFunction
	    for i = 1:length(lspace)
	    	l = lspace(i)+1;
            clf;
	    	ax = subplot(1,1,1);
	    	slice = [1:((p.nx+1) * (p.ny+1))] + l * ((p.nx+1) * (p.ny+1)); 
	    	surf(a.x, a.y, reshape(w(slice), p.ny+1, p.nx+1));
	    	xlabel('Price $p$')
	    	ylabel('Storage $q$')
            zlabel('Value')
	    	tstr = sprintf('Value for $c = %5.2f$ at $t = %4.1f$', a.z(l), p.t);
            axis([p.xmin, p.xmax, p.ymin, p.ymax])
	    	title(tstr, 'interpreter', 'latex')
	    	view([32,32])


            fname = strcat(p.outputdir,'/', p.prefix, sprintf('_value_c%d_t%d', i, p.t));
            if ~strcmp(suffix, '')
                fname = strcat(fname,'_',suffix);
            end
            set(gcf, 'Units', 'inches');
            set(gcf, 'PaperPosition', [0 0 P P]); %Position plot at left hand corner with width 5 and height 5.
            set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
            print(gcf,fname,'-dpng','-r300');
        end
    else
	    for i = 1:length(lspace)
	    	l = lspace(i)+1;
	    	ax(i) = subplot(length(lspace),1,i);
	    	slice = [1:((p.nx+1) * (p.ny+1))] + l * ((p.nx+1) * (p.ny+1)); 
	    
	    	surf(a.x, a.y, reshape(w(slice), p.ny+1, p.nx+1));
	    	if strfind(p.prefix, 'Energy')
	    		xlabel('Price p')
	    		ylabel('Storage q')
	    		tstr = sprintf('Value for $c = %4.2f$', a.z(l));
                axis([p.xmin, p.xmax, p.ymin, p.ymax])
	    		title(tstr, 'interpreter','latex')
	    		view([32,32])
	    	else
	    		tstr = sprintf('Value for $z = %4.2f$', a.z(l));
	    		title(tstr, 'interpreter','latex')
	    		xlabel('x')
	    		ylabel('y')
	    	end

	    end % for l = lspace
	    Link = linkprop(ax, {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
	    setappdata(gcf, 'StoreTheLink', Link);

    end

	

end % if p.dim == 3

end % function plotValue(p,a)
