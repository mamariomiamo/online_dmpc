% Model parameters for the (quad + controller) system
% It's model as a discrete second order system (tracking dynamics)
model_params.zeta_xy = 0.6502;
model_params.tau_xy = 0.3815;
model_params.omega_xy = 1/model_params.tau_xy;
model_params.zeta_z = 0.9103;
model_params.tau_z = 0.3;
model_params.omega_z = 1/model_params.tau_z;

% Dimension of space - 3 = 3D, 2 = 2D
ndim = 3; 

% Time settings and variables
T = 20; % simulation time 20 seconds

h = 0.2; % MPC solves the optimization every 0.2 seconds (5Hz) i.e. time step between MPC updates 0.2 seconds

tk = 0:h:20; % coarse discrete time vector for MPC updates (1*101 steps)

K = T/h+1; % number of time steps to simulate (101 steps)

Ts = 0.01; % send h/Ts = 20 commands in between MPC updates

t = 0:Ts:T; % fine simulation time vector (1*2001 steps)

% MPC related
k_hor = 16; % prediction horizon for each MPC optimization - duration of (k_hor-1)*h = 3 seconds

% Bezier-related
T_segment = 1.0; % each bezier segment has fixed duration of 1 second
deg_poly = 2; % degree of differentiability required for the position ??? why is this needed
l = 3; % number of segments of Bezier curve
d = 5; % degree of the bezier curve (5th order chosen to)

% ellipsoidal distance check
% dist_ellipsoid = E^(0.5)*dist_euclidean
c_a = [1.0, 1.0, 2.0];
ellipsoid_matrix = diag(c_a);
E1_a = ellipsoid_matrix^(-1);

% physical limits of the robot - position and acceleration bounds
phys_limits.pmin = [-1.5, -1.5, 0.2];
phys_limits.pmax = [1.5, 1.5, 2.2];
phys_limits.amax = 1;
phys_limits.amin = -1;

save('dmpc_sim_param.mat')