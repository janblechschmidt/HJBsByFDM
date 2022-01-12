
% Get derivatives of value function v, that are responsible for the optimal
% control strategy
% a python version of this file can be found in
% Dissertations/Jan\ Blechschmidt/FEniCS/pythonHJBPostprocessing
function [dec] = quadraticIntersection3d(p, c, a, v)

error('This function has to be implemented')
% 
% [~, ~, w] = determineOptimalControl(p, c, a, v, [], [], []);
% 
% dvdy = w.dvdy;
% d = 3;
% Alpha = c.alpha(a.XY);
% Beta = c.beta(a.XY);
% 
% % Initialize decision boundaries
% dec.BuyT  = zeros(0,3);
% dec.SellT = zeros(0,3);
% dec.BuyQ  = zeros(0,3);
% dec.SellQ = zeros(0,3);
% idxSell  = 1;
% 
% % Obtain triangulation of quadrilateral qrid
% tri = delaunay(a.X, a.Y, a.Z);
% numtri = size(tri,1);
% % Ensure correct dimensions
% if numtri ~= p.nx*p.ny*p.nz*6
% 	error('Wrong size of triangulation')
% end
% % triplot(tri, a.X, a.Y);
% 
% % Loop over all triangles
% for i = 1:numtri
% 	t = tri(i,:);
% 	xyt = a.XY(t,:);
% 	xt = a.X(t);
% 	yt = a.Y(t);
% 	wt = dvdy(t);
% 	at = Alpha(t);
% 	bt = Beta(t);
% 	
% 	sbuy = wt - (1 + p.costpbuy) * xt;
% 	ssell = wt - (1 - p.costpsell) * xt;
% 	csell = (ssell .* at >= p.costfsell);
% 	cbuy  = (sbuy .* bt >= p.costfbuy);
% 
% 	l = [1,1,1,2,2,3];
% 	m = [2,3,4,3,4,4];
% 
% 	if any(cbuy) ~= all(cbuy)
% 
% 		% Compute roots of quadratic polynomial along edges
% 		for k = 1:6
% 			m1 = (sbuy(l(k)) - sbuy(m(k))) * (bt(l(k)) - bt(m(k)));
% 			n = (sbuy(l(k)) - sbuy(m(k))) * bt(m(k)) + sbuy(m(k)) * (bt(l(k)) - bt(m(k)));
% 			o = sbuy(m(k)) * bt(m(k)) - p.costfbuy;
% 			if abs(m1) > 1e-8
% 				q = zeros(2,1);
% 				q(1) = (-n + sqrt(n^2 - 4*m1*o)) / (2*m1);
% 				q(2) = (-n - sqrt(n^2 - 4*m1*o)) / (2*m1);
% 			else
% 				q = - o / n;
% 			end % if abs(m1) > 1e-8
% 
% 			lidx = (q>=0) & (q<=1);
% 			if any(lidx)
% 				q0(k) = q(lidx);
% 			else
% 				q0(k) = -inf;
% 			end
% 		end
% 
% 		idx = find(q0 >= 0 & q0 <= 1);
% 
% 		switch mod(sum(cbuy),2)
% 		case 0
% 			% Quadrilateral detected
% 			% Sorting is important, otherwise matlab will plot two triangles
% 			if allOrNothing(cbuy', [1, 0, 0, 1])
% 				tmp = idx(3);
% 				idx(3) = idx(4);
% 				idx(4) = tmp;
% 			else
% 				if allOrNothing(cbuy', [1, 1, 0, 0])
% 					tmp = idx(3);
% 					idx(3) = idx(4);
% 					idx(4) = tmp;
% 				else
% 					if allOrNothing(cbuy', [1, 0, 1, 0])
% 						tmp = idx(3);
% 						idx(3) = idx(4);
% 						idx(4) = tmp;
% 					else
% 						error('no case\n')
% 					end
% 				end
% 			end
% 
% 			Q = diag(1.-q0(idx)) * xyt(l(idx),:) + diag(q0(idx)) * xyt(m(idx),:);
% 			dec.BuyQ = [dec.BuyQ; Q];
% 		case 1
% 			% Triangle detected
% 			T = diag(1.-q0(idx)) * xyt(l(idx),:) + diag(q0(idx)) * xyt(m(idx),:);
% 			dec.BuyT = [dec.BuyT; T]
% 
% 		end % switch mod(sum(cbuy),2)
% 	end % if any(cbuy) ~= all(cbuy)
% 
% 
% 	if any(csell) ~= all(csell)
% 
% 		% Compute roots of quadratic polynomial along edges
% 		for k = 1:6
% 			m1 = (ssell(l(k)) - ssell(m(k))) * (at(l(k)) - at(m(k)));
% 			n = (ssell(l(k)) - ssell(m(k))) * at(m(k)) + ssell(m(k)) * (at(l(k)) - at(m(k)));
% 			o = ssell(m(k)) * at(m(k)) - p.costfsell;
% 			if abs(m1) > 1e-8
% 				q = zeros(2,1);
% 				q(1) = (-n + sqrt(n^2 - 4*m1*o)) / (2*m1);
% 				q(2) = (-n - sqrt(n^2 - 4*m1*o)) / (2*m1);
% 			else
% 				q = - o / n;
% 			end % if abs(m1) > 1e-8
% 
% 			lidx = (q>=0) & (q<=1);
% 			if any(lidx)
% 				lidx
% 				q0(k) = q(lidx);
% 			else
% 				q0(k) = -inf;
% 			end
% 		end
% 		idx = find(q0 >= 0 & q0 <= 1)
% 		% keyboard
% 
% 		switch mod(sum(csell),2)
% 		case 0
% 			% Quadrilateral detected
% 			% Sorting is important, otherwise matlab will plot two triangles
% 			if allOrNothing(csell', [1, 0, 0, 1])
% 				tmp = idx(3);
% 				idx(3) = idx(4);
% 				idx(4) = tmp;
% 			else
% 				if allOrNothing(csell', [1, 1, 0, 0])
% 					tmp = idx(3);
% 					idx(3) = idx(4);
% 					idx(4) = tmp;
% 				else
% 					if allOrNothing(csell', [1, 0, 1, 0])
% 						tmp = idx(3);
% 						idx(3) = idx(4);
% 						idx(4) = tmp;
% 					else
% 						error('no case\n')
% 					end
% 				end
% 			end
% 
% 			Q = diag(1.-q0(idx)) * xyt(l(idx),:) + diag(q0(idx)) * xyt(m(idx),:);
% 			dec.SellQ = [dec.SellQ; Q];
% 		case 1
% 			% Triangle detected
% 			T = diag(1.-q0(idx)) * xyt(l(idx),:) + diag(q0(idx)) * xyt(m(idx),:);
% 			dec.SellT = [dec.SellT; T]
% 
% 		end % switch mod(sum(csell),2)
% 	end % if any(csell) ~= all(csell)
% 
% 
% end % for i = 1:numtri
