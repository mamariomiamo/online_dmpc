clc
clear all
close all
warning('off','all')

% Load simulation parameters
load('sim_params.mat')

% Load MPC tuning parameters
load('mpc_params.mat')

% Choose what data to visualize
global debug_constr;
debug_constr = 0;

Nplan = 10;
N = Nplan;

% for all the planning agents, extract ellipsoid-related parameters
% (Collision-related) more info about parameters see sim_params
for i = 1:Nplan
    order(i) = order_a; % order of ellipsoid
    rmin(i) = rmin_a; % X-Y protection radius for collisions
    c(i,:) = c_a;
    E1(:,:,i) = E1_a; % this is a  3*3*Nplan vector. last entry is the index of the planning agent.
    E2(:,:,i) = E2_a;
end

%% matrices used for formulation and optimization
% 1. get_model(...): model.A, model.B 
% 
% build double integrator model that represents the closed-loop dynamics of
% the position-velocity of crazyflie.
% 
% 2. get_lambda(...): Lambda (3K*3K) = (48*48)
%
% corresponds to offline-MPC equation 8: P = A0 * X0 + Lambda * U. It is
% used for applying the input and get the position sequence.
%
% 3. get_a0(...): A0 (3K*6) = (48*6) {3 axis, 16 steps, px py pz, vx, vy, vz}
% 
% Matrix represents the propagation of Initial state to the position
% sequence
%
% 4. mat_bernstein2power(...): Beta (54*54)
% 
% Dimension of Beta is tricky:
% One dimension of one segment of nth order Bezier: (n+1)*(n+1)
% 3D 3 segments 5th order: 3*3*(5+1) = 54
% convert the converts control points into polynomial coefficients. 
% P(t) = [1, t, t^2, ..., t^5] * Beta * C
% C =  [c0;c1;c2;...;c_5]
%
% 5. mat_sample_poly(...) Gamma
%
% Sample the polynomial at different time steps
%
% 6. mat_sum_sqrd_derivatives(...) Alpha
%
% Matrix used for calculate energy cost (acceleration sqr). This is the Q
% matrix in the min_snap implementation:
% min \int{(a*x)^2} = \int{x'*a'*a*x} = x' \int{a'*a} x = x'Qx
% a = [0, 0, 2, 6*t, 12*t^2, 20*t^3]
% 
% 7. H_snap: Hessian of cost function
%  H_snap = Beta'*Alpha*Beta;
% min \int{(a*Beta*C)^2} = C'Beta' \int{a'*a} BetaC = C'*Beta'*Q*Beta*C = C'*Beta'*Alpha*Beta*C
% where C is the optimizing variable a.k.a control points
%
% 8. Error to Goal cost terms
% Three senarios's weight matrix (S_free, S_obs, S_repel) 
build_all_matrices;


pmin_gen = [-1.5, -1.5, 0.2];
pmax_gen = [1.5, 1.5, 2.2];

fprintf("Planning with %i vehicles\n", N);

[po, pf] = random_test(N, pmin_gen, pmax_gen, rmin, E1, order);

use_ondemand = false;
use_repel = false;
use_stateCA = false;
use_softBVC = false;
deg_poly = 2;
run_algorithm;
p=1;
q=1;
% record the metrics for comparison
bvc_tbuild(p,q) = 1000*mean2(t_build(3:end,:));
bvc_tqp(p,q) = 1000*mean2(t_qp(3:end,:));
bvc_violated(p,q) = violated;
bvc_reachedgoal(p, q) = pass;
bvc_success(p, q) = ~violated && pass;
if pass
    bvc_totdist(p, q) = sum(sum(sqrt(diff(pos_k_i_sample(1,:,:)).^2 + ...
        diff(pos_k_i_sample(2,:,:)).^2 + ...
        diff(pos_k_i_sample(3,:,:)).^2 )));
    
    for i = 1:N
        diff_goal = pos_k_i_sample(:,:,i) - repmat(pf(:,:,i),length(t),1)';
        dist_goal = sqrt(sum(diff_goal.^2,1));
        hola = find(dist_goal >= 0.1, 1, 'last');
        if isempty(hola)
            time_index(i) = 0;
        else
            time_index(i) = hola + 1;
        end
    end
    bvc_trajtime(p, q) = max(time_index)*Ts;
    
else
    bvc_totdist(p, q) = nan;
    bvc_trajtime(p, q) = nan;
    
end
