% Get derivatives of value function v, that are responsible for the optimal
% control strategy
% a python version of this file can be found in
% Dissertations/Jan\ Blechschmidt/FEniCS/pythonHJBPostprocessing
function [dec] = postprocessMinimumArrivalTime3d(p, c, a, v)

[~, ~, w] = determineOptimalControl(p, c, a, v, [], [], []);

dvdx = w.dvdx;

% Initialize decision boundaries
dec.NegT = zeros(0,p.dim);
dec.PosT = zeros(0,p.dim);
dec.NegQ = zeros(0,p.dim);
dec.PosQ = zeros(0,p.dim);

% Obtain triangulation of quadrilateral qrid
tri = delaunay(a.X, a.Y, a.Z);
numtri = size(tri,1);
% Ensure correct dimensions
if numtri ~= p.nx*p.ny*p.nz*6
	error('Wrong size of triangulation')
end

% Loop over all triangles
for i = 1:numtri
	t = tri(i,:);
	xyt = a.XY(t,:);
	wt = dvdx(t);
	
	sp = -wt - c.beta;
	sm = -wt + c.beta;
	cp = (sp >= 0);
	cm  = (sm <= 0);
	
	l = [1,1,1,2,2,3];
	m = [2,3,4,3,4,4];
	if any(cp) ~= all(cp)

		for k = 1:6
			q(k) = (-c.beta - wt(l(k))) / (wt(m(k)) - wt(l(k))); 
		end
		bidx = (q >= 0 & q<= 1);
		idx = find(bidx);
		switch mod(sum(cp),2)
		case 0
			% Quadrilateral detected
			% Sorting is important, otherwise matlab will plot two triangles
			if allOrNothing(cp', [1, 0, 0, 1])
				tmp = idx(3);
				idx(3) = idx(4);
				idx(4) = tmp;
			else
				if allOrNothing(cp', [1, 1, 0, 0])
					tmp = idx(3);
					idx(3) = idx(4);
					idx(4) = tmp;
				else
					if allOrNothing(cp', [1, 0, 1, 0])
						tmp = idx(3);
						idx(3) = idx(4);
						idx(4) = tmp;
					else
						error('no case\n')
					end
				end
			end
			Q = diag(1.-q(idx)) * xyt(l(idx),:) + diag(q(idx)) * xyt(m(idx),:);
			dec.PosQ = [dec.PosQ;Q];
		case 1
			% Triangle detected
			T = diag(1.-q(idx)) * xyt(l(idx),:) + diag(q(idx)) * xyt(m(idx),:);
			dec.PosT = [dec.PosT; T];
		end % switch mod(sum(cp),2)

	end % if any(cp) ~= all(cp)
	if any(cm) ~= all(cm)
		for k = 1:6
			q(k) = (+c.beta - wt(l(k))) / (wt(m(k)) - wt(l(k))); 
		end
		bidx = (q >= 0 & q<= 1);
		idx = find(bidx);
		switch mod(sum(cm),2)
		case 0
			% Quadrilateral detected
			% Sorting is important, otherwise matlab will plot two triangles
			if allOrNothing(cm', [1, 0, 0, 1])
				tmp = idx(3);
				idx(3) = idx(4);
				idx(4) = tmp;
			else
				if allOrNothing(cm', [1, 1, 0, 0])
					tmp = idx(3);
					idx(3) = idx(4);
					idx(4) = tmp;
				else
					if allOrNothing(cm', [1, 0, 1, 0])
						tmp = idx(3);
						idx(3) = idx(4);
						idx(4) = tmp;
					else
						error('no case\n')
					end
				end
			end
			Q = diag(1.-q(idx)) * xyt(l(idx),:) + diag(q(idx)) * xyt(m(idx),:);
			dec.NegQ = [dec.NegQ; Q];
		case 1
			% Triangle detected
			T = diag(1.-q(idx)) * xyt(l(idx),:) + diag(q(idx)) * xyt(m(idx),:);
			dec.NegT = [dec.NegT; T];
		end % switch mod(sum(cp),2)

	end % if any(cm) ~= all(cm)
end % for i = 1:numtri
