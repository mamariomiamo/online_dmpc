%% triple integrator dynamic system
% p_dot = v
% v_dot = a
% a_dot = j

% Prediction horizon
t_hor = 4;

% time-step
syms a_0 v_0 p_0 j0 j1 j2 j3 j4;
dt = 0.2;

% 

K=5;
Ta=zeros(K);
Tv=zeros(K);
Tp=zeros(K);

for i = 1:K
    Ta(i,1:i) = ones(1,i)*dt;
end

for i = 1:K
    for j = 1:i
        Tv(i,j) = (i-j+0.5)*dt^2;
    end
end

for i = 1:K
    for j = 1:i
        Tp(i,j) = ((i-j+1)*(i-j)/2+1/6)*dt^3;
    end
end

Ba = ones(K,1)*a_0;
Bv = ones(K,1)*v_0;
Bp = ones(K,1)*p_0;

for i=1:K
    Bv(i) = Bv(i) + i*dt*a_0;
    Bp(i) = Bp(i) + i*dt*v_0 + i^2/2*a_0*dt^2;
end

Ba = ones(K,1)*a_0;
J_ = [j0;j1;j2;j3;j4];
A_ = Ta*J_+Ba
V_ = Tv*J_+Bv;
P_ = Tp*J_+Bp;
