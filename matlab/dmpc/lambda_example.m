clc
close all
clear
% Load simulation parameters
load('dmpc_sim_param.mat');
[model] = dmpc_get_model(h,model_params);
A = model.A;
B = model.B;
K = k_hor;
Lambda.pos = [];
Lambda.vel = [];
prev_row = zeros(6, 3*K); 
ncols = size(B, 2);
nrows = size(B, 1);
tic
for k = 1:K
    add_B = [zeros(nrows, ncols*(k-1)), B, zeros(nrows, ncols*(K-k))];
    new_row = A*prev_row + add_B;   
    Lambda.pos = [Lambda.pos; new_row(1:3,:)];
    %Lambda.vel = [Lambda.vel; new_row(4:6,:)];
    prev_row = new_row;   
end
toc
tic
Phi = [eye(3) zeros(3)];
basis = [];
for i=1:K
    basis = [basis, Phi*A^(K-i)*B];
end

lambda2 = [];
for i=1:K
    temp = basis(:,((K-i)*3+1):(K*3));
    new_row = [temp, zeros(3,3*(K-i))];
    lambda2 = [lambda2;new_row];
end
toc