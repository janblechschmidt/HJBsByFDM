function plotControl(p, a, u, lspace)
P=20;
SAVE = false;
if p.dim == 1
	plot(a.x, p.vec2Mat(u));
	title('Control')
	xlabel('x')
	ylabel('y')
end % if p.dim == 1

if p.dim == 2
    for i = 1:p.control_dim
        if ~p.writeControlFunction
            subplot(p.control_dim, 1, i)
        else
            f = figure('visible','off');
        end
        azsurf = surf(a.x, a.y, p.vec2Mat(u(:,i)),'EdgeAlpha',0.1,'LineWidth',0.001);
        if SAVE
        figure(1), clf,
        azsurf = surf(a.x, a.y, p.vec2Mat(u(:,i)),'EdgeColor','none','LineWidth',0.001);
        view(2)
        axis equal
        axis tight
        end
        colormap viridis;
        alpha(azsurf,0.9);
        % title(sprintf('Control %2d', i))
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
        if SAVE
        cb = colorbar('eastoutside')
        cb.Label.Interpreter='latex';
        fname=sprintf('Control%dRhoPosFD.png',i)
        set(gcf, 'Units', 'inches');
        set(gcf, 'PaperPosition', [0 0 P P]);
        set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
        print(gcf,fname,'-dpng','-r300')
        end

        % For American option case
        % axis([p.xmin,p.xmax,p.ymin,p.ymax])
        % xticks([0 20 40 60 80 100 120 140])
        % yticks([0 20 40 60 80 100 120 140])
        % view(2)

        if p.writeControlFunction
            fname = strcat(p.outputdir,'/', p.prefix, sprintf('_control_%d_theta_%02d_dofs_%05d',i, p.theta*10, length(a.dof_idxs)));
            set(gcf, 'Units', 'inches');
            set(gcf, 'PaperPosition', [0 0 P P]);
            set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
            print(gcf,fname,'-dpng','-r300'); 
            % print(gcf,fname,'-dpdf','-r300'); 
        end
    end
    if p.writeControlFunction
        fname = strcat(p.outputdir,'/', p.prefix, sprintf('_control'));
        set(gcf, 'Units', 'inches');
        set(gcf, 'PaperPosition', [0 0 P 3*P]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [P 3*P]); %Set the paper to have width 5 and height 5.
        print(gcf,fname,'-dpng','-r300'); 
        keyboard
    end

end % if p.dim == 2

if p.dim == 3

	if nargin < 5
		% lspace = round(p.nz/2);
		lspace = [round(p.nz/4), round(p.nz/2), round(p.nz*3/4)];
	end
    if p.writeControlFunction
	    for i = 1:length(lspace)
	    	l = lspace(i)+1;
            clf;
	    	ax = subplot(1,1,1);
	    	slice = [1:((p.nx+1) * (p.ny+1))] + l * ((p.nx+1) * (p.ny+1)); 
    		surf(a.x, a.y, reshape(u(slice), p.ny+1, p.nx+1));
	    	xlabel('Price $p$')
	    	ylabel('Storage $q$')
            zlabel('Control $\gamma$')
	    	tstr = sprintf('Control for $c = %5.2f$ at $t = %4.1f$', a.z(l), p.t);
            axis([p.xmin, p.xmax, p.ymin, p.ymax])
	    	title(tstr, 'interpreter', 'latex')
			% view([32,32])
            view([140,40])


            fname = strcat(p.outputdir,'/', p.prefix, sprintf('_control_c%d_t%d', i, p.t));
            set(gcf, 'Units', 'inches');
            set(gcf, 'PaperPosition', [0 0 P P]); %Position plot at left hand corner with width 5 and height 5.
            set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
            print(gcf,fname,'-dpng','-r300');
        end
    else

    	for i = 1:length(lspace)
    		l = lspace(i)+1;
    		ax(i) = subplot(length(lspace),1,i);
    		slice = [1:((p.nx+1) * (p.ny+1))] + (l-1) * ((p.nx+1) * (p.ny+1));
    		surf(a.x, a.y, reshape(u(slice), p.ny+1, p.nx+1));
    		tstr = sprintf('Control for fixed z = %4.2f', a.z(l));
    		title(tstr)
    		xlabel('x')
    		ylabel('y')
    		if strfind(p.prefix, 'Energy')
    			% hold on
    			% plot([c.meanPrice(p.t), c.meanPrice(p.t)],[p.ymin, p.ymax], 'r')
    			xlabel('Price p')
    			ylabel('Storage q')
    			tstr = sprintf('Control for fixed consumption c = %4.2f', a.z(l));
    			title(tstr)
                axis([p.xmin, p.xmax, p.ymin, p.ymax])
    			% view([32,32])
                view(2)
    
    		end
    		% view([0,90])
        end % for i = 1:length(lspace)
    end
    Link = linkprop(ax, {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
    setappdata(gcf, 'StoreTheLink', Link);

end % if p.dim == 3

end % function plotControl(p,a)
