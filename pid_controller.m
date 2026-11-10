function tau = pid_controller(mu, mu_d, psi, psi_d, params)
% ============================================================
% 修正版 PID 控制器（用于对比实验）
% 输入：mu [u;v;w;r], mu_d 期望速度, psi 当前航向, psi_d 期望航向
% 输出：tau = [tau_u; 0; tau_w; tau_r]
% ============================================================

persistent integral_psi integral_u
if isempty(integral_psi)
    integral_psi = 0; 
    integral_u = 0;
end

% 航向误差（归一化到 -pi 到 pi）
e_psi = psi_d - psi;
while e_psi > pi,  e_psi = e_psi - 2*pi; end
while e_psi < -pi, e_psi = e_psi + 2*pi; end

% 速度误差
e_u = mu_d(1) - mu(1);

% PID 参数（保守，避免震荡）
Kp_psi = 20; Ki_psi = 0.5; Kd_psi = 15;
Kp_u   = 40; Ki_u   = 1;

% 积分限幅
integral_psi = max(min(integral_psi + e_psi*params.Ts, 1), -1);
integral_u   = max(min(integral_u + e_u*params.Ts, 1), -1);

% 控制量计算（Kd 项用角速度 r 提供阻尼）
tau_r = Kp_psi*e_psi + Ki_psi*integral_psi - Kd_psi*mu(4);
tau_u = Kp_u*e_u     + Ki_u*integral_u;
tau_w = 0;

% 限幅
tau_u = max(min(tau_u, params.tau_max(1)), params.tau_min(1));
tau_r = max(min(tau_r, params.tau_max(3)), params.tau_min(3));

tau = [tau_u; 0; tau_w; tau_r];
end