function [u, ind, A, xi] = determineOptimalControl(p, a, v, vnp1)
% determineOptimalControl is a function that computes the optimal control
% for a Hamilton-Jacobian-Bellman equation
% INPUT:  p - problem structure
%         a - grid structure 
%         v - value function
%
% OUTPUT: u - optimal control for given v
% Load helper functions

helperFunctions;
m = p.ndofs;
n = p.ndofs;
EPS=1e-13;
% Overwrite system matrix with identity indicator
xi = 0;

if p.trivialControl
	u   = zeros(p.ndofs, p.control_dim);
	ind = zeros(p.ndofs, p.control_dim);
else
    % TODO: Ich lege die Differentialoperatoren aktuell erstmal hier an, eventuell lohnt es
    % sich, die einmal am Anfang bei Bedarf anzulegen

    uzero = zeros(p.ndofs, p.control_dim);
    A = sparse(p.ndofs, p.ndofs);
    for i = 1:p.control_dim

	    % Set up differential operators for computation of optimal control
	    Aforw{i} = setupDifferentialOperator(p, a, 'FirstOrderForward', i);
		Aback{i} = setupDifferentialOperator(p, a, 'FirstOrderBackward', i);

        mu = @(u) p.convection{i}(a.XY, p.t, u);
        Azero = @(u) spdiags(posPart(mu(u)),0,m,n) * Aforw{i}  - ...
            spdiags(negPart(mu(u)),0,m,n) * Aback{i};
        val_func = @(u) Azero(u) * (p.theta*v+(1-p.theta)*vnp1) + p.f(a.XY, u);
        % val_func = @(u) Azero(u) * v + p.f(a.XY, u);

        ustar_forw = zeros(p.ndofs, p.control_dim);
        ustar_back = zeros(p.ndofs, p.control_dim);
        if p.alpha > 0
            ustar_forw(:, i) = (-Aforw{i} * v - p.beta) / p.alpha;
            ustar_back(:, i) = (-Aback{i} * v + p.beta) / p.alpha;
        else
            ustar_forw(:, i) = p.Umax;
            ustar_back(:, i) = p.Umin;
        end
        % ustar_forw
        % keyboard


        % Truncate controls
        ustar_forw(ustar_forw < 0) = 0;
        ustar_back(ustar_back > 0) = 0;

        % V = [val_func(ustar_forw), val_func(ustar_back), val_func(uzero)];
        % [vmin, ind] = min(V, [], 2);
        % u(:, i) = uzero(:, i) .* (ind==3) + ustar_forw(:, i) .* (ind == 1) + ustar_back(:, i) .* (ind == 2);
        % if p.t == 0
        %     fprintf('here...')
        %     keyboard
        % end
        vzero = val_func(uzero);
        vforw = val_func(ustar_forw);
        vback = val_func(ustar_back);
        V = [vzero, vforw, vback];
        [vmin, ind] = min(V, [], 2);

        ind(abs(vmin - vforw)<EPS)=2;
        ind(abs(vmin - vback)<EPS)=3;
        ind(abs(vmin - vzero)<EPS)=1;
        u(:, i) = uzero(:, i) .* (ind==1) + ustar_forw(:, i) .* (ind == 2) + ustar_back(:, i) .* (ind == 3);

        A = A + Azero(u);
    end
end % if p.trivialControl
