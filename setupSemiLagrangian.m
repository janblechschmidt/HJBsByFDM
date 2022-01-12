function a = setupSemiLagrangian(p, a)
% setupSemiLagrangian initializes the semi-Lagrangian solver for HJB equations
% The first part can be used in all cases, the pre-computation of the matrices T and vectors F
% is only valid, if theta = 1 (fully implicit) and the time-step size is equal


if ~isfield(p, 'slMode')
    p.slMode = 1;
end

if p.trivialControl
    p.slMode = 4;
end

% slMode = 1; % only bang-bang controls alpha, 0, beta
% slMode = 2; % an equidistant grid ranging from alpha to beta (including 0, if possible)
% slMode = 3; % an equidistant grid ranging from alpha to beta (not including 0)
% slMode = 4; % trivial control, i.e. buy current consumption

% 1st step: Compute for each node i,j the best control from
switch p.slMode
case 1
    a.cntrls = p.cntrls(a.XY);
	% a.cntrls = [p.Umin, zeros(size(p.Umin)), p.Umax];
case 2
	nl = 10;
	levels = linspace(0,1,nl+1);
	% Initialize all containing u = 0
	a.cntrls = [p.Umin * (1 - levels), p.Umax * levels(2:end)];

	% Treat those for which Umin > 0 seperately
	idx = find(p.Umin > 0);
	if length(idx) > 0
		levels = linspace(0,1,2*nl+1);
		a.cntrls(idx,:) = p.Umin(idx) * (1 - levels) + p.Umax(idx) * levels;
	end

case 3
	nl = 20;
	levels = linspace(0,1,nl);
	a.cntrls = p.Umin * (1 - levels) + p.Umax * levels;
case 4
	a.cntrls = zeros(p.ndofs,1);
	% a.cntrls = a.Z;
end % switch p.slMode
n_cntrls = size(a.cntrls,2);

if p.timeConstantConvection
    switch p.dim
    case 2
    
    	[X, Y] = meshgrid(a.x, a.y);
    	myReshape = @(x) reshape(x,[p.ny+1,p.nx+1]);
    
    	for k = 1:n_cntrls
    		cntrl = myReshape(a.cntrls(:,k));
    
            t = sparse(0, 0);
    		for xi = 1:p.nx+1
    		    s = sparse(p.ny+1,p.ny+1);
    		    for yi = 1:p.ny+1
    		    	basisfunc = zeros(p.ny+1, 1);
    		    	basisfunc(yi) = 1;
    		    	tmp = interp1(Y(:,xi), basisfunc, Y(:,xi) + p.ht * p.convection{2}([X(:,xi),Y(:,xi)],p.t,cntrl(:,xi)));
    		    	% tmp = interp1(Y(:,1), basisfunc, Y(:,1) + p.ht * cntrl(:,1));
    		    	s(:,yi) = tmp(:);
    		    end % for yi = 1:p.ny+1
    			t = blkdiag(t, s);
    		end % for xi = 1:p.nx
    		a.T{k} = t;
    
    		% Get F for each part
    		a.F{k} = p.ht * p.f(a.XY, cntrl(:));
    	end % for k = 1:n_cntrls
    
    case 3
    
    	[X, Y, Z] = meshgrid(a.x, a.y, a.z);
    	myReshape = @(x) reshape(x,[p.ny+1,p.nx+1]);
    
    	for k = 1:n_cntrls
    
    		for zi = 1:(p.nz+1)
    
    			Zfix = a.z(zi);
    			Zfixidx = find(a.iz == zi);
    
    			cntrl = myReshape(a.cntrls(Zfixidx,k));
    
    			X = myReshape(a.X(Zfixidx));
    			Y = myReshape(a.Y(Zfixidx));
    			Umin = myReshape(p.Umin(Zfixidx));
    			Umax = myReshape(p.Umax(Zfixidx));
    
    			s  = sparse(p.ny+1,p.ny+1);
                fprintf('Update this')
                keyboard

    			for yi = 1:p.ny+1
    				basisfunc = zeros(p.ny+1, 1);
    				basisfunc(yi) = 1;
                    keyboard
    				tmp = interp1(Y(:,1), basisfunc, Y(:,1) + p.ht * p.convection{2}([],p.t,cntrl(:,xi)));
                    % tmp = interp1(Y(:,1), basisfunc, Y(:,1) + p.ht * cntrl(:,1) - p.ht * Zfix);
    				s(:,yi) = tmp(:);
    			end % for yi = 1:p.ny+1
    
    			t = s;
    
    			for xi = 1:p.nx
    				t  = blkdiag(t, s);
    			end % for xi = 1:p.nx
    
    			a.T{k}(Zfixidx, Zfixidx) = t;
    			a.F{k}(Zfixidx, 1) = p.ht * p.f(a.XY(Zfixidx,:),cntrl(:));
    
    		end % for zi = 1:(p.nz+1)
    
    	end % for k = 1:n_cntrls
    
    end % switch p.dim
end % p.timeConstantConvection

end % function a = setupSemiLagrangian(p, a, c)
