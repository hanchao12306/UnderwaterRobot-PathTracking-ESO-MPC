clear; clc; close all;
clear functions;
params = set_params();

% ---------- 三维螺旋路径 ----------
path = [];
for x = 0:0.5:20
    y = 1.5 * sin(0.15 * x);   % 水平面正弦
    z = 1.0 * sin(0.10 * x);   % 深度正弦
    path = [path; x, y, z];
end

eta = [0; 0; 0; 0];
mu  = [0.3; 0; 0; 0];
psi_continuous = 0; target_idx = 1; Delta = 2.5;
xi = [mu; zeros(3,1)];
f1_hat = mu;
f2_hat = zeros(4,1);
tau_prev = zeros(4,1);
rho_hat_filt = zeros(4,1);

N = 950;    % 100秒（速度0.25，路径20m）
log_eta = zeros(4, N); log_rho = zeros(4, N); log_rho_hat = zeros(4, N);
log_t = zeros(1, N); log_tau = zeros(3, N);
log_wd = zeros(1, N);   % 记录期望垂向速度

R_accept = 1.5;
tau_rho_filt = 0.2;
wd_filt = 0;            % wd 低通滤波状态
tau_w_filt = 1.0;       % wd 滤波时间常数

for k = 1:N
    tk = (k-1)*params.Ts;
    log_eta(:,k) = eta; log_t(k) = tk;
    
    % ---------- LOS 制导 ----------
    [psi_d, psi_continuous, target_idx] = ...
        los_guidance(eta, mu, path, target_idx, psi_continuous, R_accept, Delta);
    
    % ---------- 垂向制导（带低通滤波）----------
    dists_xy = sqrt((path(:,1) - eta(1)).^2 + (path(:,2) - eta(2)).^2);
    [~, nearest_idx] = min(dists_xy);
    z_ref = path(nearest_idx, 3);    % 用最近点的 z，而不是 target_idx 的 z
    wd_raw = params.kz * (z_ref - eta(3));
    wd_raw = max(min(wd_raw, params.wmax), -params.wmax);
    wd_filt = wd_filt + params.Ts/tau_w_filt * (wd_raw - wd_filt);
    wd = wd_filt;
    log_wd(k) = wd;
    
    % ---------- 纵荡速度 ----------
    ud = 0.2;
    
    % ---------- 航向误差 → r_d ----------
    e_psi = psi_d - eta(4);
    while e_psi > pi,  e_psi = e_psi - 2*pi; end
    while e_psi < -pi, e_psi = e_psi + 2*pi; end
    r_d = 0.10 * e_psi;
    r_d = max(min(r_d, 0.2), -0.2);
    
    mu_d = [ud; 0; wd; r_d];
    
    % ---------- ESO 扰动估计 ----------
    [mu_hat, rho_hat, f1_hat, f2_hat] = ...
        eso_estimator(mu, tau_prev, f1_hat, f2_hat, params);
    log_rho_hat(:,k) = rho_hat;
    
    % ---------- MPC 优化 ----------
    [tau_bar_nom, xi] = mpc_controller(mu_hat, mu_d, xi, params);
    
    % ---------- ESO 前馈补偿（不含转向通道）----------
    tau_bar_full = [tau_bar_nom(1); 0; tau_bar_nom(2); tau_bar_nom(3)];
    soft_start = min(tk / 3.0, 1.0);
    rho_hat_filt = rho_hat_filt + params.Ts/tau_rho_filt * (rho_hat - rho_hat_filt);
    rho_hat_comp = soft_start * rho_hat_filt;
    rho_hat_comp = max(min(rho_hat_comp, 0.03), -0.03);
    rho_hat_comp(3) = 0;
    rho_hat_comp(4) = 0;
    tau_full = tau_bar_full - params.M_nom * rho_hat_comp;
    tau = [tau_full(1); 0; tau_full(3); tau_full(4)];
    
    % ---------- 控制量限幅 ----------
    tau(1) = max(min(tau(1), params.tau_max(1)), params.tau_min(1));
    tau(3) = max(min(tau(3), params.tau_max(2)), params.tau_min(2));
    tau(4) = max(min(tau(4), params.tau_max(3)), params.tau_min(3));
    tau_prev = tau;
    log_tau(:,k) = tau([1,3,4]);
    
    % ---------- 真实扰动 ----------
    [~, rho_true] = auv_dynamics(mu, tau, tk, params);
    log_rho(:,k) = rho_true;
    
    % ---------- 动力学更新 ----------
    [mu, ~] = auv_dynamics(mu, tau, tk, params);
    
    % ---------- 运动学更新 ----------
    psi = eta(4);
    J = [cos(psi), -sin(psi), 0, 0;
         sin(psi),  cos(psi), 0, 0;
         0, 0, 1, 0;
         0, 0, 0, 1];
    eta = eta + params.Ts * J * mu;
end


% ----------  RMS 计算（按 X 位置匹配）----------
errors_x = zeros(1, N);
errors_y = zeros(1, N);
errors_z = zeros(1, N);
for k = 1:N
    % 找 x 坐标最接近的参考点
    [~, idx] = min(abs(path(:,1) - log_eta(1,k)));
    errors_x(k) = abs(path(idx,1) - log_eta(1,k));
    errors_y(k) = abs(path(idx,2) - log_eta(2,k));
    errors_z(k) = abs(path(idx,3) - log_eta(3,k));
end

rms_x = sqrt(mean(errors_x.^2));
rms_y = sqrt(mean(errors_y.^2));
rms_z = sqrt(mean(errors_z.^2));
rms_3d = sqrt(rms_x^2 + rms_y^2 + rms_z^2);

fprintf('==========================================\n');
fprintf('X 方向 RMS: %.4f m\n', rms_x);
fprintf('Y 方向 RMS: %.4f m\n', rms_y);
fprintf('Z 方向 RMS: %.4f m\n', rms_z);
fprintf('三维 RMS: %.4f m\n', rms_3d);
fprintf('Z 方向最大误差: %.4f m\n', max(errors_z));
fprintf('==========================================\n');

% 稳态 RMS（跳过前 10 秒）
errors_steady = sqrt(errors_x(100:end).^2 + errors_y(100:end).^2 + errors_z(100:end).^2);
fprintf('稳态 RMS（t>10s）: %.4f m\n', sqrt(mean(errors_steady.^2)));

% ---------- 三维轨迹 ----------
figure;
plot3(path(:,1), path(:,2), path(:,3), 'k--', 'LineWidth', 1.5); hold on;
plot3(log_eta(1,:), log_eta(2,:), log_eta(3,:), 'b-', 'LineWidth', 1.2);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend('参考路径', 'ESO-MPC'); grid on; view(45, 30);
title('三维路径跟踪（含海流扰动）');

% ---------- 二维投影 ----------
figure;
subplot(1,2,1);
plot(path(:,1), path(:,2), 'k--', 'LineWidth', 1.5); hold on;
plot(log_eta(1,:), log_eta(2,:), 'b-', 'LineWidth', 1.2);
xlabel('X'); ylabel('Y'); legend('参考','ESO-MPC'); grid on; axis equal;
title('水平面 XY 投影');

subplot(1,2,2);
plot(path(:,1), path(:,3), 'k--', 'LineWidth', 1.5); hold on;
plot(log_eta(1,:), log_eta(3,:), 'b-', 'LineWidth', 1.2);
xlabel('X'); ylabel('Z'); legend('参考','ESO-MPC'); grid on;
title('垂直面 XZ 投影');

% ---------- 期望速度 ----------
figure;
plot(log_t, log_wd, 'r-', 'LineWidth', 1.2);
xlabel('时间 (s)'); ylabel('w_d (m/s)');
grid on; title('期望垂向速度');