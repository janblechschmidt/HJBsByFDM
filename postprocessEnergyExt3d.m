% Get derivatives of value function v, that are responsible for the optimal
% control strategy
% a python version of this file can be found in
% Dissertations/Jan\ Blechschmidt/FEniCS/pythonHJBPostprocessing
function [dec] = postprocessEnergy3d(p, a, v)

A = setupDifferentialOperator(p, a, 'FirstOrder', 2);

dvdy = A*v;
d = 3;
Alpha = p.alpha(a.XY);
Beta = p.beta(a.XY);

% Initialize decision boundaries
dec.MinT  = zeros(0,3);
dec.MinQ = zeros(0,3);
dec.SellT = zeros(0,3);
dec.SellQ = zeros(0,3);
dec.BuyT  = zeros(0,3);
dec.BuyQ  = zeros(0,3);

% Loop over all triangles
for i = 1:p.numtri
	t = p.tri(i,:);
	xyt = a.XY(t,:);
	xt = a.X(t);
	vq = dvdy(t);
	gmin = Alpha(t);
	gmax = Beta(t);
    
	bbuy = vq - (1 + p.costpbuy) * xt - p.costfbuy;
	bsell = vq - (1 - p.costpsell) * xt + p.costfsell;

	% Boundary Do-nothing
	sect = linearIntersection3d(xyt, gmin, 0.0);
	switch size(sect,1)
	case 3
		dec.MinT = [dec.MinT; sect];
	case 4
		dec.MinQ = [dec.MinQ; sect];
	end

	% Boundary Sell
	gmin(gmin >= 0) = 0;
	sect = linearIntersection3d(xyt, bsell.*gmin, p.costcsell);
	switch size(sect,1)
	case 3
		dec.SellT = [dec.SellT; sect];
	case 4
		dec.SellQ = [dec.SellQ; sect];
	end
	
	% Boundary Buy
    if min(gmin)<0
	    sectbuy = linearIntersection3d(xyt, bbuy.*gmax, p.costcbuy);
	    switch size(sectbuy,1)
	    case 3
	    	dec.BuyT = [dec.BuyT; sectbuy];
	    case 4
	    	dec.BuyQ = [dec.BuyQ; sectbuy];
        end
    else
        % if bbuy.*gmax < p.costcbuy
    	    sectNo = linearIntersection3d(xyt, bbuy, 0.0);
    	    switch size(sectNo,1)
	        case 3
	        	dec.BuyT = [dec.BuyT; sectNo];
    	    case 4
	        	dec.BuyQ = [dec.BuyQ; sectNo];
	        end
        % end

    end

end % for i = 1:p.numtri
