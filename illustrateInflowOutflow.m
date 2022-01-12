% This script illustrates the inflow into the storage and outflow from the storage
% for the energy storage problem in 2d
clear all, clc

% Load helper functions such as positive, negative part
helperFunctions

p.dim = 2;
p.mode = 'EnergyStorage';

% Load problem specific files
p = setupParameters(p);
c = setupCoefficients(p);

N = 1000;
Q_range = linspace(p.ymin, p.ymax, N + 1)';
figure(1), clf
plot(Q_range, [c.alpha([Q_range,Q_range]), c.beta([Q_range,Q_range])] );


% Now we assume different situations
% Goal: Integrate the storage as accurate as possible, if we know, that
% the control is maximal/minimal from time t_0 to t_1

t0 = 0;
t1 = 100;
Q_t0 = 300000;


% First method: flow rate is constant
Q_t1_lin = Q_t0 + (t1-t0) * c.alpha([0,Q_t0]);
lin.t = [t0, t1];
lin.Q = [Q_t0, Q_t1_lin];

% Second method: Forward Euler (explicit)
n = 1000
traj.t = linspace(t0,t1,n+1)
traj.Q = zeros(n+1,1);
traj.Q(1) = Q_t0;
for i =1:n
	traj.Q(i+1) = traj.Q(i) + (t1-t0) / n * c.alpha([0,traj.Q(i)]);
end

% Third method: Integrate function alpha analytically
p.alpha_int = @(q0,t) ( sqrt(q0) - p.K1/2 * t ).^2 ;

figure(2), clf, hold on
plot(lin.t, lin.Q, 'r-');
plot(traj.t, traj.Q, 'k-');
plot(traj.t, p.alpha_int(Q_t0, traj.t), 'g-'); 


% Inflow...

% First method: flow rate is constant
Q_t1_lin = Q_t0 + (t1-t0) * c.beta([0,Q_t0]);
lin.t = [t0, t1];
lin.Q = [Q_t0, Q_t1_lin];

% Second method: Forward Euler (explicit)
n = 1000
traj.t = linspace(t0,t1,n+1)
traj.Q = zeros(n+1,1);
traj.Q(1) = Q_t0;
for i =1:n
	traj.Q(i+1) = traj.Q(i) + (t1-t0) / n * c.beta([0,traj.Q(i)]);
end

% Third method: Integrate function alpha analytically
% ... I found no simple way to express the integral of this function analytically

figure(3), clf, hold on
plot(lin.t, lin.Q, 'r-');
plot(traj.t, traj.Q, 'k-');

