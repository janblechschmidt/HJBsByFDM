% This script tries to simulate the value of the original optimal policy for 
% the energy storage problem and compares this to the postprocessed one.
rng(1)
% TODO: Get policy from file
if ~exist('p')
	p = setupParameters(p);
	c = setupCoefficients(p);
	[a, p] = setupStructure(p);
end
nsim = 1;
ntimesteps = p.nt;
writePlots = 0;
plotMode = 0;

% Evaluate meshgrid for interp3 command
[X, Y, Z] = meshgrid(a.x, a.y, a.z);

% Initial state
x0 = [40, 110000, 8000];
Vx0 = interp3(X, Y, Z, p.vec2Mat(V(:,1)), x0(1), x0(2), x0(3));

dt = (p.Tmax-p.Tmin) / ntimesteps;

% if ~isfield(p, 'Amin')
% 	p = setupCompareValues(p, a, c);
% end

Ay = setupDifferentialOperator(p, a, 'FirstOrder', 2);

% Maximum rates for buying and selling
Alpha = c.alpha(a.XY);
Beta = c.beta(a.XY);

% Common time index
T = a.t;

% Common mean price and consumption
muP = c.meanPrice(T);
muC = c.meanConsumption(T);

% Index matrix
% 0 - do nothing
% -1 - sell max
% +1 - buy minimum
% +2 - buy maximum
% A = zeros(size(U));
% 
% for i = 1:size(U,2)
% 	A(abs(U(:,i) - Alpha) < 1e-3 & Alpha > 0) = 1;
% 	A(abs(U(:,i) - Alpha) < 1e-3 & Alpha < 0) = -1;
% 	A(abs(U(:,i) - Beta ) < 1e-3) = 2;
% end % for i = 1:size(U,2)


for sim = 1:nsim
	if mod(sim,10) == 1
		fprintf('\n---------------------------------------------------------------------\n');
		fprintf(' Simulation |    V(x_0)   | V_orig(x_0) |  V_post(x_0) |   d_abs(x0)  | d_rel(x0)\n', sim, nsim);
		fprintf('---------------------------------------------------------------------\n');
	end

	% Common price and consumption processes
	P = zeros(1, ntimesteps);
	C = zeros(1, ntimesteps);

	% Differing processes for bought/sold energy
	S.OrigOptimal = zeros(1, ntimesteps+1);
	S.OrigPost = zeros(1, ntimesteps+1);
	S.BuyConsumption = zeros(1, ntimesteps+1);

	% Differing storage level processes
	Q.OrigOptimal = zeros(1, ntimesteps);
	Q.OrigPost = zeros(1, ntimesteps);
	Q.BuyConsumption = zeros(1, ntimesteps);

	% Differing cost/gain processes
	F.OrigPost       = zeros(1, ntimesteps);
	F.OrigOptimal    = zeros(1, ntimesteps);
	F.BuyConsumption = zeros(1, ntimesteps);

	% Simulate paths of Brownian motions for price and consumption processes
	xi = randn(2, ntimesteps);

	% Simulate path of price process for this simulation
	P(1) = x0(1);
	for i = 1:ntimesteps
		P(i+1) = P(i) + p.kappap * ( muP(i) - P(i)) * dt + p.sigmap*sqrt(dt) * xi(1,i);
	end % for i = 1:ntimesteps

	% Simulate path of consumption process for this simulation
	C(1) = x0(3);
	for i = 1:ntimesteps
		C(i+1) = C(i) + p.kappac * ( muC(i) - C(i)) * dt + sqrt(dt) * p.sigmac *(p.correlationcp * xi(1,i) + sqrt(1-p.correlationcp^2) * xi(2,i));
	end % for i = 1:ntimesteps

	% Initialize paths for different policies for storage level process
	Q.OrigOptimal(1) = x0(2);
	Q.OrigPost(1) = x0(2);
	Q.BuyConsumption(1) = x0(2);

	% ---------------------------------------------
	% 1st version: Buy consumption (trivial policy)
	% ---------------------------------------------
	for i = 1:ntimesteps

		S.BuyConsumption(i) = C(i);
		F.BuyConsumption(i) = c.f([P(i), Q.BuyConsumption(i), C(i)], S.BuyConsumption(i)) * dt;
		% Q.BuyConsumption(i+1) = Q.BuyConsumption(i) + (S.BuyConsumption(i) - C(i))*dt;
		Q.BuyConsumption(i+1) = Q.BuyConsumption(i);

	end % for i = 1:ntimesteps
	F.BuyConsumption(ntimesteps + 1) = c.finalTimeVal([P(i+1), Q.BuyConsumption(i+1), C(i+1)]);


	% ---------------------------------------------
	% 2nd version: Original optimal policy (not postprocessed)
	% ---------------------------------------------
	for i = 1:ntimesteps

		% Make sure that the storage  level doesn't become negative
		if Q.OrigOptimal(i) < 0
			if Q.OrigOptimal(i) > - 1e-8
				Q.OrigOptimal(i) = 0;
			else
				keyboard
			end % if Q.OrigOptimal(i) > - 1e-8
		end % if Q.OrigOptimal(i) < 0

		% We are given the current Price, Consumption and Storage Level
		% and compute the current policy
		S.OrigOptimal(i) = interp3(X, Y, Z, p.vec2Mat(U(:,i)), P(i), Q.OrigOptimal(i), C(i));

		F.OrigOptimal(i)    = c.f([P(i), Q.OrigOptimal(i),    C(i)], S.OrigOptimal(i)) * dt;

		Q.OrigOptimal(i+1) = Q.OrigOptimal(i) + (S.OrigOptimal(i) - C(i)) * dt;

	end % for i = 1:ntimesteps
	F.OrigOptimal(ntimesteps + 1) = c.finalTimeVal([P(i+1), Q.OrigOptimal(i+1), C(i+1)]);

	% -----------------------------------------
	% 3rd version: Postprocessed optimal policy
	% -----------------------------------------
	for i = 1:ntimesteps

		if Q.OrigPost(i) < 0
			if Q.OrigPost(i) > - 1e-8
				Q.OrigPost(i) = 0;
			else
				keyboard
			end % if Q.OrigPost(i) > - 1e-8
		end % if Q.OrigPost(i) < 0

		mode = 'FOO'; % First-order optimality
		% mode = 'CV';  % Compare values
		switch mode
		case 'FOO'
			wt = interp3(X, Y, Z, p.vec2Mat(Ay*V(:,i)), P(i), Q.OrigPost(i), C(i));
			
			at = max( [c.alpha([P(i), Q.OrigPost(i), C(i)]) , p.qmin + C(i) - Q.OrigPost(i) ]);
			bt = min( [c.beta([P(i), Q.OrigPost(i), C(i)])  , p.qmax + C(i) - Q.OrigPost(i) ]);

			ct = at * (at > 0); % Minimal available control


			
			ssell = wt - (1 - p.costpsell) * P(i);
			sbuy = wt - (1 + p.costpbuy) * P(i);
			
			sellopt = (ssell * at >= p.costfsell);
			buyopt = (sbuy * bt >= p.costfbuy);
			if sellopt
				ct = at;
			end
			if buyopt
				ct = bt;
			end

			if i == ntimesteps
				tmp = c.fillLevel + C(i) - Q.OrigPost(i);
				if at <= tmp & tmp <= bt 
					ct = tmp;
				end
				% if ct - C(i) + Q.OrigPost(i) < c.fillLevel
				% 	tmp = c.fillLevel + C(i) - Q.OrigPost(i);
				% 
				% 	if at <= tmp & tmp <= bt 
				% 		ct = tmp;
				% 	end
				% end
			end

		case 'CV'

			M = cell2mat(cellfun(@(x,y) x*V(:,i)+y , a.T, a.F,'UniformOutput',false));
			[rhs, idx] = max(M,[],2);
			u = a.cntrls(sub2ind(size(a.cntrls),(1:size(idx,1))', idx ));

			keyboard
			% Vumin  = p.fmin + p.Amin*V(:,i);
			% Vumax  = p.fmax + p.Amax*V(:,i);
			% Vuzero = p.Azero*V(:,i);
			% wsell = interp3(X,Y,Z,p.vec2Mat(Vumin), P(i), Q.OrigPost(i), C(i));
			% wbuy  = interp3(X,Y,Z,p.vec2Mat(Vumax), P(i), Q.OrigPost(i), C(i));
			% wzero = interp3(X,Y,Z,p.vec2Mat(Vuzero), P(i), Q.OrigPost(i), C(i));
			[wmin, wind] = max([wsell, wzero, wbuy],[],2);
			switch wind
			case 1
				ct = at;
			case 2
				ct = 0;
			case 3
				ct = bt;
			end % switch wind
		end % switch mode
	
		S.OrigPost(i) = ct;

		F.OrigPost(i)       = c.f([P(i), Q.OrigPost(i),       C(i)], S.OrigPost(i)) * dt;
	
		% Update storage level for different versions
		Q.OrigPost(i+1) = Q.OrigPost(i) +  (S.OrigPost(i) - C(i)) * dt;


	end % for i = 1:ntimesteps
	F.OrigPost(ntimesteps + 1) = c.finalTimeVal([P(i+1), Q.OrigPost(i+1), C(i+1)]);

		% Postprocessed
		% tri = p.Tri.pointLocation(P(i),Q.OrigOptimal(i),C(i))
		% vertices = p.Tri.ConnectivityList(tri,:)
		% xyt = a.XY(vertices,:);
		% xt = a.X(vertices);
		% at = Alpha(vertices);
		% bt = Beta(vertices);
		% if i == 143
		% 	keyboard
		% end

	
		% if Q.OrigPost(i) - C(i) + at < p.qmin
		% 	at = C(i) - Q.OrigPost(i);
		% end % if Q.OrigPost(i) - C(i) + at < p.qmin
		% 
		% if Q.OrigPost(i) - C(i) + bt > p.qmax
		% 	bt = p.qmax + C(i) - Q.OrigPost(i);
		% end % if Q.OrigPost(i) - C(i) + bt > p.qmax

		% Vumin  = p.fmin + p.Amin*V(:,i);
		% Vumax  = p.fmax + p.Amax*V(:,i);
		% Vuzero = p.Azero*V(:,i);

		% tmp = [Vumin, Vumax, Vuzero];
		% [~,wind] = max(tmp,[],2);

		% wzero = interp3(X,Y,Z,p.vec2Mat(V(:,i)), P(i), Q.OrigPost(i) - C(i), C(i));
		% wsell = interp3(X,Y,Z,p.vec2Mat(V(:,i)), P(i), Q.OrigPost(i) - C(i) + at, C(i));
		% wbuy  = interp3(X,Y,Z,p.vec2Mat(V(:,i)), P(i), Q.OrigPost(i) - C(i) + bt, C(i));
		% wsell = interp3(X,Y,Z,p.vec2Mat(Vumin), P(i), Q.OrigPost(i) - C(i) + at, C(i));
		% wbuy  = interp3(X,Y,Z,p.vec2Mat(Vumax), P(i), Q.OrigPost(i) - C(i) + bt, C(i));
		% wzero = interp3(X,Y,Z,p.vec2Mat(Vuzero), P(i), Q.OrigPost(i) - C(i) + 0, C(i));
	
		% 
		% ct = at * (at > 0); % Minimal available control
		% 
		% ssell = wt - (1 - p.costpsell) * P(i);
		% sbuy = wt - (1 + p.costpbuy) * P(i);
		% 
		% sellopt = (ssell * at >= p.costfsell);
		% buyopt = (sbuy * bt >= p.costfbuy);
		% if sellopt
		% 	ct = at;
		% end
		% if buyopt
		% 	ct = bt;
		% end
		% % ct = at * (sellopt) + bt * (buyopt);
		% if Q.OrigPost(i) +  (S.OrigPost(i) - C(i))*dt < 0
		% 	S.OrigPost(i) = -Q.OrigPost(i)/(dt) +C(i);
		% else
		% 	S.OrigPost(i) = ct;
		% end
	
	% keyboard
	F_BuyConsumption(sim) = sum(F.BuyConsumption);
	F_Post(sim) = sum(F.OrigPost);
	F_Orig(sim) = sum(F.OrigOptimal);
	d_abs(sim) = F_Post(sim) - F_Orig(sim);
	d_rel(sim) = d_abs(sim) / abs(F_Orig(sim));

	fprintf('%4d / %4d | %6.4e | %6.4e | %6.4e  |  %6.4e  |  %6.4f \n', sim, nsim, Vx0, F_Orig(sim), F_Post(sim), d_abs(sim), d_rel(sim));
end

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
	mean(F_BuyConsumption)
end

keyboard
% Store simulated data as csv
DATA = [a.t', P', C', Q.BuyConsumption', Q.OrigOptimal', Q.OrigPost', S.BuyConsumption', S.OrigOptimal', S.OrigPost',  F.BuyConsumption', F.OrigOptimal', F.OrigPost']  
csvwrite([p.outputdir '/' p.prefix 'simulatedData.csv'], DATA);

if plotMode == 1

	sfigure(4); clf

	subplot(3,1,1)
	plot(T,P,'r-'); hold on
	plot(T,c.meanPrice(T),'r-','LineWidth',0.5)
	axis([p.Tmin, p.Tmax, p.pmin, p.pmax] )
	title('Price process')

	subplot(3,1,2)
	plot(T,Q.OrigPost,'mo-'), hold on
	plot(T,Q.OrigOptimal,'b+-') 
	plot(T,Q.BuyConsumption,'r+-')
	axis([p.Tmin, p.Tmax, p.qmin, p.qmax] )
	title('Storage process')

	subplot(3,1,3)
	plot(T,C,'b-'); hold on
	plot(T,c.meanConsumption(T),'b-','LineWidth',0.5)
	axis([p.Tmin, p.Tmax, p.cmin, p.cmax] )
	title('Consumption process')

	sfigure(5); clf;
	plot(T(1:end-1),S.OrigOptimal,'b-'), hold on
	plot(T(1:end-1),S.OrigPost,'m-')
	axis([p.Tmin, p.Tmax, -20000, 20000] )
	title('Controlled process')
end % if plotMode == 1

if plotMode == 2
	figure(1); clf

	plot(T,P,'r-'); hold on
	plot(T,c.meanPrice(T),'r-','LineWidth',0.5)
	axis([p.Tmin, p.Tmax, p.pmin, p.pmax] )
	title('Price process')
	xlabel('Time')
	if writePlots
		print([p.outputdir '/example_PriceProcess'], '-dpdf')
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
		print([p.outputdir '/example_StorageProcess'], '-dpdf')
	end

	figure(3); clf
	plot(T,C,'k-'); hold on
	plot(T,c.meanConsumption(T),'k-','LineWidth',0.5)
	axis([p.Tmin, p.Tmax, p.cmin, p.cmax] )
	title('Consumption process')
	xlabel('Time')
	if writePlots
		print([p.outputdir '/example_ConsumptionProcess'], '-dpdf')
	end

	figure(4); clf
	plot(T(1:end-1),S.OrigPost,'m-'), hold on
	plot(T(1:end-1),S.OrigOptimal,'b-')
	axis([p.Tmin, p.Tmax, -20000, 20000] )
	title('Controlled process')
	xlabel('Time')
	legend('Postprocessed', 'Original', 'Location', 'southwest')
	if writePlots
		print([p.outputdir '/example_OptimalControl'], '-dpdf')
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
			print(sprintf([p.outputdir '/control_%d_%03d'],p.nz, i),'-dpng')
		end
	end
end % if plotMode == 3
