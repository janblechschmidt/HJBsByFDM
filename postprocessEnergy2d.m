% Get derivatives of value function v, that are responsible for the optimal
% control strategy
% a python version of this file can be found in
% Dissertations/Jan\ Blechschmidt/FEniCS/pythonHJBPostprocessing
function [dec] = postprocessEnergy2d(p, c, a, v)

[~, ~, w] = determineOptimalControl(p, c, a, v, [], [], []);

dvdy = w.dvdy;

Alpha = c.alpha(a.XY);
Beta = c.beta(a.XY);

% Initialize decision boundaries
dec.Buy  = zeros(0,p.dim);
dec.Sell = zeros(0,p.dim);
idxBuy  = 1;
idxSell  = 1;

% Loop over all triangles
for i = 1:p.numtri
	t = p.tri(i,:);
	xyt = a.XY(t,:);
	xt = a.X(t);
	yt = a.Y(t);
	wt = dvdy(t);
	at = Alpha(t);
	bt = Beta(t);
	
	sbuy = wt - (1 + p.costpbuy) * xt;
	ssell = wt - (1 - p.costpsell) * xt;
	csell = (ssell .* at >= p.costfsell);
	cbuy  = (sbuy .* bt >= p.costfbuy);

	if any(cbuy) ~= all(cbuy)
		for k = 1:2
			for l = (k+1):3
				if cbuy(l) ~= cbuy(k)
					m = (sbuy(l) - sbuy(k)) * (bt(l) - bt(k));
					n = (sbuy(l) - sbuy(k)) * bt(k) + sbuy(k) * (bt(l) - bt(k));
					o = sbuy(k) * bt(k) - p.costfbuy;
					if abs(m) > 1e-8
						q = zeros(2,1);
						q(1) = (-n + sqrt(n^2 - 4*m*o)) / (2*m);
						q(2) = (-n - sqrt(n^2 - 4*m*o)) / (2*m);
					else
						q = - o / n;
					end % if abs(m) > 1e-8

					q0 = q((q>=0) & (q<=1));

					if length(q0) ~= 1
						error('something with computation of root is wrong')
					end

					% Append this point to decision boundary
					dec.Buy(idxBuy,:) = (1.-q0) * xyt(k,:) + q0 * xyt(l,:);
					idxBuy = idxBuy + 1;

				end % if cbuy(l) != cbuy(k):
			end % for l = (k+1):3
		end % for k = 1:2
	end % if any(cbuy) ~= all(cbuy)

	if any(csell) ~= all(csell)
		for k = 1:2
			for l = (k+1):3
				if csell(l) ~= csell(k)
					m = (ssell(l) - ssell(k)) * (at(l) - at(k));
					n = (ssell(l) - ssell(k)) * at(k) + ssell(k) * (at(l) - at(k));
					o = ssell(k) * at(k) - p.costfsell;
					if abs(m) > 1e-8
						q = zeros(2,1);
						q(1) = (-n + sqrt(n^2 - 4*m*o)) / (2*m);
						q(2) = (-n - sqrt(n^2 - 4*m*o)) / (2*m);
					else
						q = - o / n;
					end % if abs(m) > 1e-8

					q0 = q((q>=0) & (q<=1));

					if length(q0) ~= 1
						error('something with computation of root is wrong')
					end

					% Append this point to decision boundary
					dec.Sell(idxSell,:) = (1.-q0) * xyt(k,:) + q0 * xyt(l,:);
					idxSell = idxSell + 1;

				end % if csell(l) != csell(k):
			end % for l = (k+1):3
		end % for k = 1:2
	end % if any(csell) ~= all(csell)

end % for i = 1:p.numtri
% 
% % Plot decision boundary
% sfigure(3), clf, hold on
% plotDecisionBoundary(dec.Buy, 'k-');
% plotDecisionBoundary(dec.Sell, 'r-');
