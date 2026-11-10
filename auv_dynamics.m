function [mu_next, rho_true] = auv_dynamics(mu, tau, t, params)
% ============================================================
% AUV 4自由度动力学模型
% 对应论文式(2-6)：M*mu_dot + C(mu)*mu + D(mu)*mu + g(eta) = tau + tau_d
% 输入：mu  = [u;v;w;r] 当前速度
%       tau = [tau_u; 0; tau_w; tau_r] 控制输入（4维，tau_v=0）
% 输出：mu_next 下一时刻速度
%       rho_true 真实集总扰动（用于与ESO估计对比）
% ============================================================

M    = params.M;
Minv = params.Minv;
C    = params.C_func(mu);    % 科氏力矩阵（式2-9）
D    = params.D_func(mu);    % 阻尼矩阵（式2-11）
g    = zeros(4,1);           % 恢复力，假设中性浮力（假设2.3）

% 外部扰动力/力矩（式2-6 中的 tau_d）
tau_d = params.tau_d_func(t);

% 海流影响（等效为附加扰动力）
Vc = params.current_func(t);
tau_current = [Vc(1); Vc(2); Vc(3); 0];

% ---------- 真实动力学（式2-6） ----------
% mu_dot = M^{-1} * (tau + tau_d + tau_current - C*mu - D*mu - g)
mu_dot  = Minv * (tau + tau_d + tau_current - C*mu - D*mu - g);
mu_next = mu + params.Ts * mu_dot;

% ---------- 真实集总扰动 rho（式2-17） ----------
% rho = M^{-1}*(-C*mu - D*mu - g + tau_d + tau_current) + (M^{-1} - M_nom^{-1})*tau
rho_true = Minv * (-C*mu - D*mu - g + tau_d + tau_current) ...
         + (Minv - params.Minv_nom) * tau;

end