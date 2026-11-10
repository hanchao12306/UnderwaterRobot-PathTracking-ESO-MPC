function [tau_bar_nom, xi] = mpc_controller(mu, mu_d, xi_prev, params)
% ============================================================
% MPC 控制器（论文4.4节）
% 名义系统（式4-35）：mu_dot = -M_nom^{-1}*D_nom*mu + M_nom^{-1}*E43*tau_bar
% 增广状态空间（式4-38）：xi = [mu; tau_bar]
% ============================================================

Ts = params.Ts; Np = params.Np; Nc = params.Nc;
Q  = params.Q;  R  = params.R;

% E43 矩阵：将3维控制量映射到4维速度通道（v通道无控制）
E43 = [1 0 0;
       0 0 0;
       0 1 0;
       0 0 1];

% 名义阻尼矩阵
D_nom = diag([params.d11, params.d22, params.d33, params.d44]);

% 【关键修正】正确的离散化公式，必须包含 M_nom^{-1}
A_nom = eye(4) - Ts * (params.Minv_nom * D_nom); 
B_nom = Ts * params.Minv_nom * E43;   

% 增广状态空间模型（式4-38）
A_aug = [A_nom, B_nom;
         zeros(3,4), eye(3)];
B_aug = [B_nom;
         eye(3)];
C_aug = [eye(4), zeros(4,3)];

% 预测方程（式4-44）
Psi   = zeros(4*Np, 7);
Theta = zeros(4*Np, 3*Nc);
A_pow = eye(7);
for i = 1:Np
    A_pow = A_pow * A_aug;
    Psi((i-1)*4+1:i*4, :) = C_aug * A_pow;
    for j = 1:min(i,Nc)
        A_pow2 = eye(7);
        for kk = 1:i-j
            A_pow2 = A_pow2 * A_aug;
        end
        Theta((i-1)*4+1:i*4, (j-1)*3+1:j*3) = C_aug * A_pow2 * B_aug;
    end
end

Y_ref = repmat(mu_d, Np, 1);
Q_bar = kron(eye(Np), Q);
R_bar = kron(eye(Nc), R);

xi_prev(1:4) = mu;

H = Theta' * Q_bar * Theta + R_bar;
H = (H + H')/2;
f = Theta' * Q_bar * (Psi * xi_prev - Y_ref);

% 控制增量约束（式4-48）
du_min = repmat(params.du_min, Nc, 1);
du_max = repmat(params.du_max, Nc, 1);
A_ineq = [ eye(3*Nc); -eye(3*Nc)];
b_ineq = [ du_max;   -du_min];

options = optimoptions('quadprog','Display','off');
[dU, ~, exitflag] = quadprog(H, f, A_ineq, b_ineq, [], [], [], [], [], options);
if exitflag ~= 1, dU = zeros(3*Nc, 1); end

% 取第一个控制增量（式4-65）
dU_first     = dU(1:3);
tau_bar_prev = xi_prev(5:7);
tau_bar_nom  = tau_bar_prev + dU_first;

xi = xi_prev;
xi(5:7) = tau_bar_nom;
end