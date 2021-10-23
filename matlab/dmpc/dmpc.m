clc;
close all;
clear;

%% simulation parameters (omitted VICON noise and some of the ellipsoid params)
% Load simulation parameters
load('dmpc_sim_param.mat');
% Load MPC optimization parameters
load('dmpc_mpc_param.mat');

%% MPC
% 1. continuous model -> discrete model 
% (in this case is the p-v tracking dynamics parameterized by model.A and
% model.B.
% 
% 2. discrete model -> prediction model
% position sequence P{} for 16 prediction steps for 3 axis
% Get prediction matrix, 

%% buid matrices

% discrete model
% x_{n+1} = A * x_n + B * u_n
% x is system state [x y z vx vy vz]'
% u is position reference [x y z]'
% model.A is 6*6
% model.B is 6*3
[model] = dmpc_get_model(h, model_params);

% Lambda
% corresponds to offline-MPC equation 8: P = A0 * X0 + Lambda * U. It is
% used for applying the input and get the position sequence.
% Lambda = [ Phi*B,         zeros(3), ..., zeros(3);
%                      Phi*A*B,     Phi*B,       ]
Lambda = dmpc_get_lambda(model.A, model.B, k_hor);

% A0 reflects the propagatino of initial state
A0 = dmpc_get_a0(model.A, k_hor);

% Beta - converts control points into polynomial coefficients
Beta = mat_bernstein2power(d, l, ndim);

% Gamma - sample the polynomial at different time steps
Gamma = mat_sample_poly(T_segment, 0:h:((k_hor-1)*h), d, l);

% Alpha - sum of squared derivatives of the position
cr = zeros(1, d+1); % weights on degree of derivative deg = [0 1 ... d]
cr(3) = cost_acc;
Alpha = mat_sum_sqrd_derivatives(d, T_segment, cr, l, ndim);

%%%%%%%%%%%%%% CONSTRUCT TERMS OF THE COST FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Hessian for the minimum snap cost function is
H_snap = Beta'*Alpha*Beta;

% For the goal tracking error cost function, define a weight matrix S
% Case 1: no collisions in the horizon - go fast
S_free = s_free*[zeros(3*(k_hor-spd_f), 3*k_hor);
                 zeros(3*spd_f, 3*(k_hor-spd_f)) eye(3*spd_f)];
             
% Case 2: collisions in the horizon - go slower
S_obs = s_obs*[zeros(3*(k_hor-spd_o), 3*k_hor);
               zeros(3*spd_o, 3*(k_hor-spd_o)) eye(3*spd_o)];

% Case 3: colliding scenario, change setpoint and repel rapidly
S_repel = s_repel*[eye(3*spd_r) zeros(3*spd_r, 3*(k_hor-spd_r));
                   zeros(3*(k_hor-spd_r), 3*k_hor)];
               
Phi = Lambda.pos*Gamma*Beta; %% ??? how is this derived (eq12 from offline DMPC) convert bezier to natural basis then sample
Phi_vel = Lambda.vel*Gamma*Beta; %% Gamma is sampling from Bezier, Beta is BezierToPolynomial
H_free = Phi'*S_free*Phi;
H_obs = Phi'*S_obs*Phi;
H_repel = Phi'*S_repel*Phi;

Phi_ref = Gamma*Beta;

% The complete Hessian is simply the sum of the two
H_f = H_free + H_snap;
H_o = H_obs + H_snap;
H_r = H_repel + H_snap;

% Predefine this matrix to construct f_pf_repel when needed
Rho_repel = S_repel*Phi;

% We can also construct the matrix that will then be multiplied by the
% initial condition -> X0'*A0'*S*Lambda*Gamma*Beta
% f is the gradient term to be passed into QP solver
mat_f_x0_free = A0.pos'*S_free*Lambda.pos*Gamma*Beta;
mat_f_x0_obs = A0.pos'*S_obs*Lambda.pos*Gamma*Beta;
mat_f_x0_repel = A0.pos'*S_repel*Lambda.pos*Gamma*Beta;

%%%%%%%%%%%%% CONSTRUCT INEQUALITY CONSTRAINT MATRICES %%%%%%%%%%%%%%%%%%%%%%

% Types of constraints: 1) Acceleration limits 2) workspace boundaries
% We need to create the matrices that map position control points into c-th
% derivative control points
T_ctrl_pts = cell_derivate_ctrl_pts(d);
[A_in, b_in] = build_ineq_constraints(d, l, h, ndim, k_hor, T_segment,...
                                      phys_limits, T_ctrl_pts, Beta);

%%%%%%%%%%%%% CONSTRUCT EQUALITY CONSTRAINT MATRICES %%%%%%%%%%%%%%%%%%%%

% Types of constraints: 1) Continuity constraints up to degree deg_poly
A_eq = build_eq_constraints(d, l, ndim, deg_poly, T_ctrl_pts);
% The constant vector b_eq will be updated within the solver function,
% since it requires knowledge of the initial condition of the reference

%%%%%%%%%%%%%% MATRICES TO DECODE SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First construct all the matrices that map the solution vector to samples
% of the n-th derivative of position
[model_s, inv_model_s] = get_model(Ts, model_params);
K_sample = length(0:Ts:h-Ts);
Lambda_s = get_lambda(model_s.A, model_s.B, K_sample);
A0_s = get_a0(model_s.A, K_sample);

for r = 0:d
    if r > 0
        Mu = T_ctrl_pts{r};
        Mu_3d = augment_array_ndim(Mu, 3);
        Sigma_r = kron(eye(l), Mu_3d);
    else
        Sigma_r = eye(3*(d+1)*l);
    end
    
    Beta_r = mat_bernstein2power(d-r, l, ndim);     
    
    % Sample Bezier curves at 1/h Hz for the whole horizon
    t_sample_r = 0:h:((k_hor-1)*h);
    Tau_r = mat_sample_poly(T_segment, t_sample_r, d-r, l);
    Der_h{r+1} = Tau_r*Beta_r*Sigma_r;
    
    % Sample Bezier curves at 1/Ts Hz for the first applied input
    t_sample_r = 0:Ts:h-Ts;
    Tau_r = mat_sample_poly(T_segment, t_sample_r, d-r, l);
    Der_ts{r+1} = Tau_r*Beta_r*Sigma_r;
end

% Sample states at 1/Ts Hz for the first applied input
t_sample_r = 0:Ts:h-Ts;
Tau_r = mat_sample_poly(T_segment, t_sample_r, d, l);
Phi_sample = Lambda_s.pos*Tau_r*Beta;
Phi_vel_sample = Lambda_s.vel*Tau_r*Beta;

% workspace boundary
pmin_gen = [-1.5, -1.5, 0.2];
pmax_gen = [1.5, 1.5, 2.2];

ctr = [0;0;1]; % Center of Boundary Sphere
rad = [0.75;0.75;1]; % Radius of Boundary Sphere
zlim = 0.3; % z limit cannot be more than 1

% number of planning agents
N = 5;
rmin = 0.3; 
E1 = E1_a;
order = order_a;
dist = 1;

for i = 1:N
    order(i) = order_a;
    rmin(i) = rmin_a;
    c(i,:) = c_a;
    E1(:,:,i) = E1_a;
    E2(:,:,i) = E2_a;
end
% Generate a random set of initial and final positions
[po, pf] = random_test(N, pmin_gen, pmax_gen, rmin, E1, order);

%% Init position
for i = 1:N
    StartPose(i,:) = RandSphereCoord(rad,ctr,zlim);
    
    if i > 1
        while norm(StartPose(i,:) - StartPose(i-1,:)) < dist
            StartPose(i,:) = RandSphereCoord(rad,ctr,zlim);
        end
    end
    EndPose(i,:) = [ctr(1)-(StartPose(i,1) - ctr(1)), ...
        ctr(2)-(StartPose(i,2) - ctr(2)), ...
        ctr(3)-(StartPose(i,3) - ctr(3))]
    plot3(StartPose(i,1),StartPose(i,2),StartPose(i,3),'-o','Color','b','MarkerSize',10,...
    'MarkerFaceColor','#D9FFFF');
    hold on
    plot3(EndPose(i,1),EndPose(i,2),EndPose(i,3),'-o','Color','b','MarkerSize',10,...
    'MarkerFaceColor','#D9FFFF');
    line([StartPose(i,1),StartPose(i,2),StartPose(i,3)],[EndPose(i,1),EndPose(i,2),EndPose(i,3)],'Color', 'w');

end

%% Algo
% Run algorithm with Collision avoidance in the state space.
use_ondemand = true;
use_repel = false;
use_stateCA = true;
use_softBVC = false;
deg_poly = 2;

for i = 1:N
   f_pf_free(:,:,i) = repmat((pf(:,:,i))', k_hor, 1)'*S_free*Phi;
   f_pf_obs(:,:,i) = repmat((pf(:,:,i))', k_hor, 1)'*S_obs*Phi; 
   poi = po(:,:,i)';
   voi = 0.001*ones(3, 1);
   X0(:,i) = [poi; voi];
   pos_k_i(:,1,i) = poi;
   vel_k_i(:,1,i) = voi;
   pos_k_i_sample(:,1,i) = poi;
   X0_ref(:,:,i) = [poi, voi, zeros(3,d - 1)];
   prev_state(:,i) = X0(:,i);
   for r = 1:deg_poly+1
      ref(:,1,r,i) = X0_ref(:,r,i); 
      ref_sample(:,1,r,i) = X0_ref(:,r,i);
   end
   hor_ref(:,:,i,1) = repmat(poi, 1, k_hor);
   hor_rob(:,:,i,1) = repmat(poi, 1, k_hor+1);
end

pred_X0 = X0;

% Variables for reference replanning based on state feedback % ???
integ_err(:,1) = zeros(3, 1);

num_repels = zeros(N,1);

for k = 2:K % K is the number of MPC update steps
    for i=1:N
        % Compare the expected and sensed position at time k
        err_pos(:,k) = X0(1:3,i) - pred_X0(1:3,i);
        err_vel(:,k) = X0(4:6,i) - pred_X0(4:6,i);
        
        % Compare the current position and the reference
        err_pos_ref(:,k) = X0(1:3,i) - X0_ref(:,1,i);
        err_vel_ref(:,k) = X0(4:6,i) - X0_ref(:,2,i);
        
        der_err_ref(:,k) = (err_pos_ref(:,k) - err_pos_ref(:,k-1)) / h;
        
        % Cost that determines whether there's something disturbing the agent
        % Cost gets higher when the error between the reference and the state gets higher
        cost(:,k,i) = (err_pos_ref(:,k).^5) ./ (-X0(4:6,i)+sign(X0(4:6,i))*0.01);
        
        % Integral term on position
        integ_err(:,k) = integ_err(:,k-1) + err_pos_ref(:,k)*h;
        
        for n = 1:ndim
            % Filter noise on position for feedback
            if abs(err_pos(n,k)) < err_tol_pos
                err_pos(n,k) = 0;
            end
        end
        
        trigger(k,i) = 0;
        
        % Reset reference to state if the error grows large
        if any(cost(:,k,i) > max_cost) || any(cost(:,k,i) < min_cost)
            X0_ref(:,1,i) = X0(1:3,i);
            X0_ref(:,2,i) = X0(4:6,i);
            X0_ref(:,3:5,i) = zeros(3,3);
            trigger(k,i) = 1;
        else
            X0_ref(:,1,i) = X0_ref(:,1,i); %+ err_pos(:,k) + ki*integ_err(:,k);
        end
        
        % Include on-demand collision avoidance
        [A_coll, b_coll, pf_tmp, t_build(k,i)] = ondemand_softconstraints(hor_rob(:,2:end,:,k-1), Phi,...
            X0(:,i), A0.pos, i, rmin,...
            order, E1, E2);
        if ~isempty(b_coll) % collisions in the horizon
            % Include collision constraints and slack variables
            N_v = length(b_coll) / 3;
            A_in_i = [A_in zeros(size(A_in,1), N_v) ; A_coll];
            b_in_i = [b_in; b_coll];
            A_eq_i = [A_eq zeros(size(A_eq,1), N_v)];
            
            % Linear and quadratic term to penalize collision relaxation
            f_eps = lin_coll_penalty*ones(1, N_v);
            H_eps = quad_coll_penalty*eye(N_v);
            
            % If close to colliding, change setpoint to quickly react
            if ~isempty(pf_tmp) && use_repel
                num_repels(i) = num_repels(i) + 1;
                H_i = [H_r zeros(size(H_f,1), N_v);
                    zeros(N_v, size(H_f,2)) H_eps];
                mat_f_x0_i = mat_f_x0_repel;
                f_tot = repmat((pf_tmp),k_hor,1)'*Rho_repel;
            else
                H_i = [H_o zeros(size(H_f,1), N_v);
                    zeros(N_v,size(H_f,2)) H_eps];
                mat_f_x0_i = mat_f_x0_obs;
                f_tot = f_pf_obs(:,:,i);
            end
            
        else % no collisions in horizon
            A_in_i = A_in;
            b_in_i = b_in;
            A_eq_i = A_eq;
            H_i = H_f;
            f_eps = [];
            mat_f_x0_i = mat_f_x0_free;
            f_tot = f_pf_free(:,:,i);
        end
        % Solve QP
        t_start = tic;
        [sol, exitflag] = softMPC_update(l, deg_poly, A_in_i, b_in_i, A_eq_i, H_i,...
            mat_f_x0_i, f_tot, f_eps, X0(:,i), X0_ref(:,:,i));
        
        t_qp(k,i) = toc(t_start);
        
        
        %% Post-processing
        if  isempty(sol)
            x = prev_x{i};
            %             assert(~isempty(x), 'ERROR: No solution found - exitflag =  %i\n',exitflag);
        else
            prev_x{i} = sol;
            x = sol;
        end
        
        %% 
        
        
        
        % Extract the control points
        u = x(1:size(mat_f_x0_free, 2));
     
        % Apply input to model starting form our previous init condition
        pos_i = vec2mat(Phi*u + A0.pos*X0(:,i),3)';
        vel_i = vec2mat(Phi_vel*u + A0.vel*X0(:,i),3)';
        
        % Sample at a higher frequency the interval 0:Ts:h-Ts
        % This tells us what should be the value of our state after
        % sending the optimal commands if the model was perfect
        pos_i_sample = vec2mat(Phi_sample*u + A0_s.pos*X0(:,i),3)';
        vel_i_sample = vec2mat(Phi_vel_sample*u + A0_s.vel*X0(:,i),3)';
        
        % Sample the resulting reference Bezier curves at 1/h and 1/Ts
        % Get the next input to be applied 'X0_ref'
        cols = 2 + (k-2)*(h/Ts):1 + (k-1)*(h/Ts);
        for r = 1:d+1
            rth_ref(:,:,r) = vec2mat(Der_h{r}*u, 3)';
            rth_ref_sample(:,:,r) = vec2mat(Der_ts{r}*u, 3)';
            X0_ref(:,r,i) = rth_ref(:,2,r);
            ref(:,k,r,i) = rth_ref(:,2,r);
            ref_sample(:,cols,r,i) = rth_ref_sample(:,:,r);
        end
        
        % Simulate sending trajectories every Ts and applying at each time
        % step noise to the measurements and propagating the state forward
        X0_ex(:,1) = X0(:,i);
        for k_ex = 2:length(t_sample_r) + 1
            X0_ex(:, k_ex -1) = X0_ex(:, k_ex -1) + rnd_noise(std_p,std_v);
            X0_ex(:,k_ex) = model_s.A*X0_ex(:, k_ex-1) + model_s.B*rth_ref_sample(:, k_ex-1, 1);
        end
        
        % Initial conditions for next MPC cycle - based on sensing
        X0(:,i) = X0_ex(:, end);

        % Update agent's states at 1/h and 1/Ts frequencies
        pos_k_i_sample(:,cols,i) = X0_ex(1:3, 2:end);
        vel_k_i_sample(:,cols,i) = X0_ex(4:6, 2:end);
        
        pred_X0(:,i) = [pos_i_sample(:,end); vel_i_sample(:,end)];
        pos_k_i(:,k,i) = X0(1:3,i);
        vel_k_i(:,k,i) = X0(4:6,i);            
        
        % Reference and state prediction horizons - visualization purposes
        hor_ref(:,:,i,k) = rth_ref(:,:,1);
        hor_rob_k(:,:,i) = [X0(1:3,i) pos_i(:,1:end)];
        
        if i==0
            hor_rob_k(:,:,i) = [X0(1:3,i) repmat(X0(1:3,i),1,k_hor)];
        else
            hor_rob_k(:,:,i) = [X0(1:3,i) pos_i(:,1:end)];
        end
    end
    hor_rob(:,:,:,k) = hor_rob_k;
end
