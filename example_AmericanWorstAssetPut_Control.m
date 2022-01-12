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

% Overwrite system matrix with identity indicator
xi = 1;

if p.trivialControl
	u   = zeros(p.ndofs, 1);
	ind = zeros(p.ndofs, 1);
else

    uzero = zeros(p.ndofs, p.control_dim);
    A = a.A1func(p.t,uzero);

    % resA = v - p.boundaryVal(p.t, a.XY);
    resA = v - p.boundaryVal(p.Tmax, a.XY);
    resI = -(vnp1 - v)/p.ht;
    resI = resI - (a.A2 + A - a.A0) * ((1-p.theta)*vnp1 + p.theta*v) - p.f(a.XY, uzero);
    u = zeros(p.ndofs, 1);
    u(resA < resI) = 1;
    ind = u;
end % if p.trivialControl
