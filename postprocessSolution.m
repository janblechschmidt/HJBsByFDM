% Get derivatives of value function v, that are responsible for the optimal
% control strategy. A python version of this file can be found in
% Dissertations/Jan\ Blechschmidt/FEniCS/pythonHJBPostprocessing
function [decBndrs, h] = postprocessSolution(p, a, v, plotFlag)

if nargin < 4
	plotFlag = 0;
end

switch p.dim
case 2
	if strfind(p.prefix, 'Energy')
		decBndrs = postprocessEnergy2d(p, a, v)
		% Plot decision boundary
		if plotFlag
			h = sfigure(3); clf, hold on
			plotDecisionBoundary2d(decBndrs.Buy, 'k-');
			plotDecisionBoundary2d(decBndrs.Sell, 'r-');
		end
	end
	if strfind(p.prefix, 'MAT')
		decBndrs = postprocessMinimumArrivalTime2d(p, a, v);
		% Plot decision boundary
		if plotFlag
			h = sfigure(3); clf, hold on
			plotDecisionBoundary2d(decBndrs.Neg, 'k-');
			plotDecisionBoundary2d(decBndrs.Pos, 'r-');
		end
	end
	axis([p.xmin,p.xmax,p.ymin,p.ymax]);

case 3
	if strfind(p.prefix, 'Energy')
		decBndrs = postprocessEnergy3d(p, a, v);
		% Plot decision boundary
		if plotFlag
			h = sfigure(3); clf, hold on
			% h = figure('Visible', 'off'); clf, hold on
            % plotDecisionBoundary3d(decBndrs.MinT, 'triangle', 'k','FaceAlpha', 0.1, 'EdgeAlpha',0.05,'EdgeColor','k');
            % plotDecisionBoundary3d(decBndrs.MinQ, 'quadrilateral', 'k','FaceAlpha', 0.1, 'EdgeAlpha',0.05,'EdgeColor','k');
			plotDecisionBoundary3d(decBndrs.BuyT, 'triangle', 'g','FaceAlpha', 0.7, 'EdgeAlpha',0.5,'EdgeColor','g');
			% alpha(0.7);
			plotDecisionBoundary3d(decBndrs.BuyQ, 'quadrilateral', 'g','FaceAlpha', 0.7, 'EdgeAlpha',0.5,'EdgeColor','g');
			% alpha(0.7);
			plotDecisionBoundary3d(decBndrs.SellT, 'triangle', 'r','FaceAlpha', 0.7, 'EdgeAlpha',0.5,'EdgeColor','r');
			% alpha(0.7);
			plotDecisionBoundary3d(decBndrs.SellQ, 'quadrilateral', 'r','FaceAlpha', 0.7 ,'EdgeAlpha',0.5,'EdgeColor','r');
			% alpha(0.7);
            xticks([p.xmin, 0.5*(p.xmin + p.xmax), p.xmax])
            yticks([p.ymin, 0.5*(p.ymin + p.ymax), p.ymax])
            yticklabels([0, 0.5, 1.0])
            zticks([p.zmin, 0.5*(p.zmin + p.zmax), p.zmax])
			% xlabel('Price');
			% ylabel('Storage');
			% zlabel('Consumption');
			xlabel('$p$');
            ax = gca;
            ax.YAxis.Exponent = 0;
			ylabel('$q/q_{\mathrm{max}}$');
			zlabel('$c$');
			view([20,40])
		end
	end
	if strfind(p.prefix, 'MAT')
		decBndrs = postprocessMinimumArrivalTime3d(p, a, v);
		% Plot decision boundary
		if plotFlag
			h = sfigure(3); clf, hold on
			plotDecisionBoundary3d(decBndrs.NegT, 'triangle', 'k');
			plotDecisionBoundary3d(decBndrs.PosT, 'triangle', 'r');
			plotDecisionBoundary3d(decBndrs.NegQ, 'quadrilateral', 'k');
			plotDecisionBoundary3d(decBndrs.PosQ, 'quadrilateral', 'r');
		end
	end
	axis([p.xmin,p.xmax,p.ymin,p.ymax,p.zmin,p.zmax]);
end
% 
% % Plot decision boundary
% sfigure(3), clf, hold on
% plotDecisionBoundary(decBuy, 'k-');
% plotDecisionBoundary(decSell, 'r-');
