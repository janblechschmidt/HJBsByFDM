% Get derivatives of value function v, that are responsible for the optimal
% control strategy
% a python version of this file can be found in
% Dissertations/Jan\ Blechschmidt/FEniCS/pythonHJBPostprocessing
function [dec] = postprocessMinimumArrivalTime2d(p, c, a, v)

[~, ~, w] = determineOptimalControl(p, c, a, v, [], [], []);

dvdx = w.dvdx;

% Initialize decision boundaries
dec.Neg = zeros(0,p.dim);
dec.Pos = zeros(0,p.dim);
idxNeg  = 1;
idxPos = 1;

% Obtain triangulation of quadrilateral qrid
tri = delaunay(a.X, a.Y);
numtri = size(tri,1);
% Ensure correct dimensions
if numtri ~= p.nx*p.ny*2
	error('Wrong size of triangulation')
end
% triplot(tri, a.X, a.Y);

% Loop over all triangles
for i = 1:numtri
	t = tri(i,:);
	xyt = a.XY(t,:);
	xt = a.X(t);
	yt = a.Y(t);
	wt = dvdx(t);
	
	sp = -wt - c.beta;
	sm = -wt + c.beta;
	cp = (sp >= 0);
	cm  = (sm <= 0);

	if any(cp) ~= all(cp)

		l0 = (-c.beta - wt(1)) / (wt(2) - wt(1));
		l1 = (-c.beta - wt(1)) / (wt(3) - wt(1));
		l2 = (-c.beta - wt(2)) / (wt(3) - wt(2));
		if l0 >= 0 && l0 <= 1
			dec.Pos(idxPos,:) = (1.-l0) * xyt(1,:) + l0 * xyt(2,:);
			idxPos = idxPos + 1;
		end
		if l1 >= 0 && l1 <= 1
			dec.Pos(idxPos,:) = (1.-l1) * xyt(1,:) + l1 * xyt(3,:);
			idxPos = idxPos + 1;
		end
		if l2 >= 0 && l2 <= 1
			dec.Pos(idxPos,:) = (1.-l2) * xyt(2,:) + l2 * xyt(3,:);
			idxPos = idxPos + 1;
		end
	end % if any(cp) ~= all(cp)

	if any(cm) ~= all(cm)

		l0 = (+c.beta - wt(1)) / (wt(2) - wt(1));
		l1 = (+c.beta - wt(1)) / (wt(3) - wt(1));
		l2 = (+c.beta - wt(2)) / (wt(3) - wt(2));
		if l0 >= 0 && l0 <= 1
			dec.Neg(idxNeg,:) = (1.-l0) * xyt(1,:) + l0 * xyt(2,:);
			idxNeg = idxNeg + 1;
		end
		if l1 >= 0 && l1 <= 1
			dec.Neg(idxNeg,:) = (1.-l1) * xyt(1,:) + l1 * xyt(3,:);
			idxNeg = idxNeg + 1;
		end
		if l2 >= 0 && l2 <= 1
			dec.Neg(idxNeg,:) = (1.-l2) * xyt(2,:) + l2 * xyt(3,:);
			idxNeg = idxNeg + 1;
		end
	end % if any(cm) ~= all(cm)
end % for i = 1:numtri
