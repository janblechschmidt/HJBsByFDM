% This script tries to simulate the value of the original optimal policy for 
% the energy storage problem and compares this to the postprocessed one.

rng(1)
% if ~exist('p')
% 	p = setupParameters(p);
% 	c = setupCoefficients(p);
% 	[a, p] = setupStructure(p);
% end
% 
nsim = 100000;
ntimesteps = p.nt;
writePlots = 0;
plotMode = 0;
P = 20;

% Parameters for compare values
% Set control to postprocessed control value
mode = 'FOO'; % First-order optimality
% mode = 'CV';  % Compare values
nl = 10;
levels = linspace(0,1,nl+1);
levels0 = linspace(0,1,2*nl+1);

% Evaluate meshgrid for interp3 command
[X, Y, Z] = meshgrid(a.x, a.y, a.z);

% Initial state
x0 = [40, 110000, 8000];
Vx0 = interp3(X, Y, Z, p.vec2Mat(V(:,1)), x0(1), x0(2), x0(3));

dt = (p.Tmax-p.Tmin) / ntimesteps;

Ay = setupDifferentialOperator(p, a, 'FirstOrder', 2);

% Common time index
T = a.t;

% Common price and consumption processes
PP = zeros(nsim, ntimesteps);

Orig.Q = zeros(nsim, ntimesteps);
Post.Q = zeros(nsim, ntimesteps);
Zero.Q = zeros(nsim, ntimesteps);

CC = zeros(nsim, ntimesteps);

Orig.S = zeros(nsim, ntimesteps);
Post.S = zeros(nsim, ntimesteps);
Zero.S = zeros(nsim, ntimesteps);

Orig.F = zeros(nsim, ntimesteps);
Post.F = zeros(nsim, ntimesteps);
Zero.F = zeros(nsim, ntimesteps);

% Simulate paths of Brownian motions for price and consumption processes
eta_p = randn(nsim, ntimesteps);
eta_c = randn(nsim, ntimesteps);

% Initialize paths
PP(:,1) = x0(1);
CC(:,1) = x0(3);

Zero.Q(:,1) = x0(2);
Post.Q(:,1) = x0(2);
Orig.Q(:,1) = x0(2);


% Simulate processes for trivial control

fprintf('\n\nDetermine simulated values...\n')
fprintf('Simulate trivial control\n')
for i=1:ntimesteps
	t = a.t(i);
	xi = [PP(:,i), Zero.Q(:,i), CC(:,i)];

	% Set control to consumption
	Zero.S(:,i) = CC(:,i);

	Zero.F(:,i) = p.f(xi, Zero.S(:,i)) * dt;

	% Compute convection
	B = cellfun(@(c) c(xi, t, Zero.S(:,i)), p.convection, 'UniformOutput', false);

	PP(:,i+1) = PP(:,i) + B{1} * dt + p.sigmap*sqrt(dt) * eta_p(:,i);
	Zero.Q(:,i+1) = Zero.Q(:,i) + B{2} * dt;
	CC(:,i+1) = CC(:,i) + B{3} * dt + p.sigmac*sqrt(dt) *(p.correlationcp * eta_p(:,i) + sqrt(1-p.correlationcp^2) * eta_c(:,i));
end % for i = 1:ntimesteps
Zero.F(:,ntimesteps + 1) = p.finalTimeVal([PP(:,i+1), Zero.Q(:,i+1), CC(:,i+1)]);

% F_Zero can be used to test whether the value is indeed the real option value.
% Current output should be storage cost * x0(2) * 
F_Zero = sum(Zero.F, 2);
fprintf('mean(F_Zero): %f\n', mean(F_Zero));
fprintf('-x0(2)*p.coststorage*365: %f\n', -x0(2)*p.coststorage*365);

% Simulate original control (i.e. not postprocessed)
fprintf('Simulate original control\n')
for i = 1:ntimesteps
	t = a.t(i);
	xi = [PP(:,i), Orig.Q(:,i), CC(:,i)];
	% Set control to interpolated control value
	Orig.S(:,i) = interp3(X, Y, Z, p.vec2Mat(U(:,i)), xi(:,1), xi(:,2), xi(:,3));

	Orig.F(:,i) = p.f(xi, Orig.S(:,i)) * dt;

	% Compute convection
	B = cellfun(@(c) c(xi, t, Orig.S(:,i)), p.convection, 'UniformOutput', false);

	Orig.Q(:,i+1) = Orig.Q(:,i) + B{2} * dt;

	almost_negative_Q = find(Orig.Q(:,i+1) < 0  & Orig.Q(:,i+1) > -1e-10);
    if any(almost_negative_Q)
        % fprintf('found almost negative Q')
        % keyboard
    end
    if any(Orig.Q(:,i+1) <= -1e-10)
        fprintf('found negative Q')
        keyboard
    end
	Orig.Q(almost_negative_Q,i+1) = 0;

end % for i = 1:ntimesteps
Orig.F(:,ntimesteps + 1) = p.finalTimeVal([PP(:,i+1), Orig.Q(:,i+1), CC(:,i+1)]);
F_Orig = sum(Orig.F, 2);
fprintf('mean(F_Orig): %f\n', mean(F_Orig));

% Simulate postprocessed control
fprintf('Simulate postprocessed control\n')
for i = 1:ntimesteps
	t = a.t(i);
	xiprev = xi;
	xi = [PP(:,i), Post.Q(:,i), CC(:,i)];

	switch mode
	case 'FOO'
		wt = interp3(X, Y, Z, p.vec2Mat(Ay*V(:,i)), xi(:,1), xi(:,2), xi(:,3));

        betaplus = wt - (1+p.costpbuy) * xi(:,1) - p.costfbuy;
        betaminus = wt - (1-p.costpsell) * xi(:,1) + p.costfsell;
        % Maximum buy and sell rates might be changed to ensure discrete admissibility
        % Q_{t+1} = Q_t + dt * (gamma - C_t) \in [qmin,qmax]
        gbuy = min([p.beta(xi), CC(:,i)+(p.qmax-Post.Q(:,i))/dt],[],2);
        gsell = max([p.alpha(xi),  CC(:,i)+(p.qmin-Post.Q(:,i))/dt],[],2);

        buyopt = (gsell > 0 & betaplus >= 0) | (gsell <= 0 & gbuy .* betaplus >= p.costcbuy) ;
        sellopt =  (gsell > 0 & betaminus < 0) | (gsell <= 0 & gsell .* betaminus >= p.costcsell) ;

        ct = zeros(nsim,1);
        ct(gsell > 0) = gsell(gsell > 0);
        ct(buyopt) = gbuy(buyopt);
        ct(sellopt) = gsell(sellopt);


		% at = max( [p.alpha(xi) , p.qmin + CC(:,i) - Post.Q(:,i) ],[],2);
		% bt = min( [p.beta(xi)  , p.qmax + CC(:,i) - Post.Q(:,i) ],[],2);
		% ct = at .* (at > 0); % Minimal available control

		% % HERE I should use an inline function, the problem is that
		% % at can be positive that would mean buying, but instead I use always sell costs
		% ssell = wt - (1 - p.costpsell) * PP(:,i);
		% sbuy = wt - (1 + p.costpbuy) * PP(:,i);
		
		% buyopt = (sbuy .* bt >= p.costfbuy);
		% sellopt = (ssell .* at .* (at < 0) + sbuy .* at .* (at > 0) >= p.costfsell) & (~buyopt);

		% ct(sellopt) = at(sellopt);
		% ct(buyopt) = bt(buyopt);

		if i == ntimesteps
			tmp = CC(:,i) + (p.fillLevel - Post.Q(:,i))/dt;
			tmpIsPossible = (gsell <= tmp & tmp <= gbuy);
			ct(tmpIsPossible) = tmp(tmpIsPossible);
		end

	case 'CV'
		at = max( [p.alpha(xi) , p.qmin + CC(:,i) - Post.Q(:,i) ],[],2);
		bt = min( [p.beta(xi)  , p.qmax + CC(:,i) - Post.Q(:,i) ],[],2);
		% Initialize all containing u = 0
		cntrls = [at * (1 - levels), bt * levels(2:end)];
		
		% Treat those for which Umin > 0 seperately
		idx = find(at > 0);
		if length(idx) > 0
			% keyboard
			cntrls(idx,:) = at(idx,:) * (1 - levels0) + bt(idx,:) * levels0;
		end
		M = zeros(nsim,2*nl+1);
		for j = 1:(2*nl+1)
			M(:,j) = dt * p.f(xi, cntrls(:,j)) + interp3(X,Y,Z, p.vec2Mat(V(:,i+1)), xi(:,1), xi(:,2) + dt * (cntrls(:,j) - xi(:,3)), xi(:,3));
		end % for j = 1:(2*nl+1)
		[rhs, idx] = max(M,[],2);
		ct = cntrls(sub2ind(size(cntrls),(1:size(idx,1))', idx ));

		% at = max( [p.alpha(xi) , p.qmin + CC(:,i) - Post.Q(:,i) ],[],2);
		% bt = min( [p.beta(xi)  , p.qmax + CC(:,i) - Post.Q(:,i) ],[],2);
		% ct = at .* (at > 0); % Minimal available control
		% 
		% fval_at = p.f(xi, at);
		% fval_bt = p.f(xi, bt);
		% fval_ct = p.f(xi, ct);
		% Vval_at = interp3(X,Y,Z, p.vec2Mat(V(:,i+1)), xi(:,1), xi(:,2) + dt * (at - xi(:,3)), xi(:,3));
		% Vval_bt = interp3(X,Y,Z, p.vec2Mat(V(:,i+1)), xi(:,1), xi(:,2) + dt * (bt - xi(:,3)), xi(:,3));
		% Vval_ct = interp3(X,Y,Z, p.vec2Mat(V(:,i+1)), xi(:,1), xi(:,2) + dt * (ct - xi(:,3)), xi(:,3));
		% M = [Vval_at + fval_at, Vval_bt + fval_bt, Vval_ct + fval_ct];

		% [rhs, idx] = max(M,[],2);
		% ct = at.*(idx==1) + bt.*(idx==2) + ct.*(idx==3);
	end % switch mode
	
	Post.S(:,i) = ct;

	Post.F(:,i) = p.f(xi, ct) * dt;

	% Compute convection
	B = cellfun(@(c) c(xi, t, Post.S(:,i)), p.convection, 'UniformOutput', false);
    B2 = p.convection{2}(xi, t, ct);
	Post.Q(:,i+1) = Post.Q(:,i) + B2 * dt;
    % if any(Post.Q(:,i+1) < 0)
    %     keyboard
    % end

	almost_negative_Q = find(Post.Q(:,i+1) < 0  & Post.Q(:,i+1) > -1e-10);
    if any(almost_negative_Q)
        Post.Q(almost_negative_Q,i+1) = 0;
        % fprintf('Found almost negative Q')
        % keyboard
    end
    if any(Post.Q(:,i+1) <= -1e-10)
        fprintf('found negative Q')
        keyboard
    end

end % for i = 1:ntimesteps
Post.F(:,ntimesteps + 1) = p.finalTimeVal([PP(:,i+1), Post.Q(:,i+1), CC(:,i+1)]);

F_Post = sum(Post.F, 2);
fprintf('mean(F_Post): %f\n', mean(F_Post));
valorig = mean(F_Orig);
valpost = mean(F_Post);
d_abs = F_Post - F_Orig;
d_rel = d_abs ./ abs(F_Orig);

	% F_BuyConsumption(sim) = sum(F.BuyConsumption);
	% F_Post(sim) = sum(F.OrigPost);
	% F_Orig(sim) = sum(F.OrigOptimal);
	% d_abs(sim) = F_Post(sim) - F_Orig(sim);
	% d_rel(sim) = d_abs(sim) / abs(F_Orig(sim));

	% fprintf('%4d / %4d | %6.4e | %6.4e | %6.4e  |  %6.4e  |  %6.4f \n', sim, nsim, Vx0, F_Orig(sim), F_Post(sim), d_abs(sim), d_rel(sim));
if nsim > 5
	imp = F_Post./F_Orig - 1;
	imp(isnan(imp)) = [];
	if plotMode > 0
		figure(6), clf,
		hist(imp,50)
	end
	mean(imp)
	median(imp)
	% interp3(X, Y, Z, p.vec2Mat(V(:,1)), x0(1), x0(2), x0(3))
	Vx0
	mean(F_Orig)-Vx0
	mean(F_Post)-Vx0
	% This is equal to x0(2) * 365 * p.coststorage
	mean(F_Zero)
end

idx = 1;
if plotMode == 1

	sfigure(4); clf

	subplot(3,1,1)
	plot(T,PP(idx,:),'r-'); hold on
	plot(T,p.meanPrice(T),'r--','LineWidth',0.5)
	axis([p.Tmin, p.Tmax, p.pmin, p.pmax] )
	title('Price process')

	subplot(3,1,2)
	plot(T,Post.Q(idx,:),'m-'), hold on
	plot(T,Orig.Q(idx,:),'b-') 
	plot(T,Zero.Q(idx,:),'r-')
	axis([p.Tmin, p.Tmax, p.qmin, p.qmax] )
	title('Storage process')

	subplot(3,1,3)
	plot(T, CC(idx,:),'b-'); hold on
	plot(T,p.meanConsumption(T),'b--','LineWidth',0.5)
	axis([p.Tmin, p.Tmax, p.cmin, p.cmax] )
	title('Consumption process')

	sfigure(5); clf;
	plot(T(1:end-1),Orig.S(idx,:),'b-'), hold on
	plot(T(1:end-1),Post.S(idx,:),'m-')
	axis([p.Tmin, p.Tmax, -30000, 30000] )
	title('Controlled process')
end % if plotMode == 1

if plotMode == 2
	figure(1); clf

	plot(T,P,'r-'); hold on
	plot(T,p.meanPrice(T),'r-','LineWidth',0.5)
	axis([p.Tmin, p.Tmax, p.pmin, p.pmax] )
	title('Price process')
	xlabel('Time')
	if writePlots
        set(gcf, 'Units', 'inches');
        set(gcf, 'PaperPosition', [0 0 P P]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
        fname = sprintf('%s/example_PriceProcess', p.outputdir);
        print(gcf, fname,'-dpng','-r300'); 
		% print('', '-dpdf')
	end

	figure(2); clf
	plot(T,Q.OrigPost,'mo-'), hold on
	plot(T,Q.OrigOptimal,'b+-') 
	% plot(T,Q.BuyConsumption,'r+-')
	axis([p.Tmin, p.Tmax, p.qmin, p.qmax] )
	title('Storage process')
	xlabel('Time')
	legend('Postprocessed', 'Original', 'Location', 'southwest')
	if writePlots
        set(gcf, 'Units', 'inches');
        set(gcf, 'PaperPosition', [0 0 P P]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
        fname = sprintf('%s/example_StorageProcess', p.outputdir);
        print(gcf, fname,'-dpng','-r300'); 
	end

	figure(3); clf
	plot(T,C,'k-'); hold on
	plot(T,p.meanConsumption(T),'k-','LineWidth',0.5)
	axis([p.Tmin, p.Tmax, p.cmin, p.cmax] )
	title('Consumption process')
	xlabel('Time')
	if writePlots
        set(gcf, 'Units', 'inches');
        set(gcf, 'PaperPosition', [0 0 P P]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
        fname = sprintf('%s/example_ConsumptionProcess', p.outputdir);
        print(gcf, fname,'-dpng','-r300'); 
	end

	figure(4); clf
	plot(T(1:end-1),S.OrigPost,'m-'), hold on
	plot(T(1:end-1),S.OrigOptimal,'b-')
	axis([p.Tmin, p.Tmax, -20000, 20000] )
	title('Controlled process')
	xlabel('Time')
	legend('Postprocessed', 'Original', 'Location', 'southwest')
	if writePlots
        set(gcf, 'Units', 'inches');
        set(gcf, 'PaperPosition', [0 0 P P]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [P P]); %Set the paper to have width 5 and height 5.
        fname = sprintf('%s/example_OptimalControl', p.outputdir);
        print(gcf, fname,'-dpng','-r300'); 
	end
end % if plotMode == 2

if plotMode == 3
	% Plot unpostprocessed control
	for i = 1:p.nz+1
		figure(1), clf,
		plotControl(p, a, u, i)
		% view([160,20])
		view([0,90])
		xlabel('Price')
		ylabel('Storage')
		caxis([-20000, 20000])
		title(sprintf('Control for fixed consumption c = %6.2f',a.z(i)))
		if writePlots
			print(sprintf('%s/control_%d_%03d', p.outputdir, p.nz, i),'-dpng')
		end
	end
end % if plotMode == 3
