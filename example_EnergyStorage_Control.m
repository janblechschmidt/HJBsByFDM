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
xi = 0;

Aqforw = setupDifferentialOperator(p, a, 'FirstOrderForward', 2);
Aqback = setupDifferentialOperator(p, a, 'FirstOrderBackward', 2);

conv = @(u) p.convection{2}(a.XY, p.t, u);
Aq = @(u) spdiags(posPart(conv(u)),0,m,n) * Aqforw  - ...
			spdiags(negPart(conv(u)),0,m,n) * Aqback;

val_func = @(u) Aq(u) * ((1-p.theta)*vnp1 + p.theta*v) + p.f1(a.XY, u);
% val_func = @(u) Aq(u) * v + p.f1(a.XY, u);
Umin = p.alpha(a.XY);
Umax = p.beta(a.XY);
Uzero = zeros(size(Umin));
Uzero(Umin>0) = Umin(Umin>0);

V = [val_func(Umin), val_func(Umax), val_func(Uzero)];
[vmax, ind] = max(V, [], 2);
u = Umin .* (ind==1) + Umax .* (ind==2) + Uzero .* (ind==3);
A = a.A1func(p.t, u);

end
