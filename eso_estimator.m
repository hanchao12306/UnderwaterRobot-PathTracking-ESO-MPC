function [mu_hat, rho_hat, f1_hat, f2_hat] = ...
    eso_estimator(mu, tau, f1_hat, f2_hat, params)
% ============================================================
% 线性扩张状态观测器 LESO（论文4.3节）
% 系统: mu_dot = M_nom^{-1}*tau + rho（式2-27）
% 扩张状态: f1 = mu, f2 = rho（4维扩张到8维）
% 增益: l1i = 2*wo, l2i = wo^2（式4-25 带宽参数化）
% ============================================================

Ts = params.Ts; 
wo = params.omega_o;

L1 = 2*wo * eye(4);
L2 = wo^2 * eye(4);
Minv_nom = params.Minv_nom;

% 观测误差
e = mu - f1_hat;

% 前向欧拉离散化（式4-24 的离散实现）
f1_hat_next = f1_hat + Ts * (Minv_nom * tau + f2_hat + L1 * e);
f2_hat_next = f2_hat + Ts * (L2 * e);

mu_hat  = f1_hat_next;
rho_hat = f2_hat_next;
f1_hat  = f1_hat_next;
f2_hat  = f2_hat_next;
end