clear; clc; close all;
clear functions;
utils_plot();
params = set_params();

% ============================================================
% 加载 A*+APF 规划路径
% ============================================================
if exist('astar_path.mat', 'file')
    load('astar_path.mat', 'path', 'start_pos', 'goal_pos');
    fprintf('已加载 A*+APF 规划路径，点数: %d\n', size(path,1));
    fprintf('起点: (%.2f, %.2f, %.2f)\n', start_pos(1), start_pos(2), start_pos(3));
    fprintf('终点: (%.2f, %.2f, %.2f)\n', goal_pos(1), goal_pos(2), goal_pos(3));
else
    error('未找到 astar_path.mat，请先运行规划代码并保存');
end

% 强制保证路径点数 ≥ 200
if size(path, 1) < 200
    t = 1:size(path,1);
    tq = linspace(1, size(path,1), 300);
    path = [interp1(t, path(:,1), tq, 'spline')', ...
            interp1(t, path(:,2), tq, 'spline')', ...
            interp1(t, path(:,3), tq, 'spline')'];
end

% 计算路径总长度，估算仿真时长
path_len = 0;
for i = 1:size(path,1)-1
    path_len = path_len + norm(path(i+1,:) - path(i,:));
end
fprintf('路径总长度: %.2f m\n', path_len);

ud = 0.1;                                   % 巡航速度 m/s
T_est = path_len / ud + 20;                  % 估算时长 + 20s 冗余
N = ceil(T_est / params.Ts);                 % 仿真步数
fprintf('预计仿真时长: %.1f s，步数: %d\n', T_est, N);

% ============================================================
% AUV 初始化（从路径起点出发）
% ============================================================
eta = [path(1,1); path(1,2); path(1,3); 0];

% 初始航向：取路径前 5 个点的方向
psi0 = atan2(path(5,2)-path(1,2), path(5,1)-path(1,1));
eta(4) = psi0;

mu  = [0.3; 0; 0; 0];
psi_continuous = psi0;
target_idx = 1;
Delta = 0.8;

xi = [mu; zeros(3,1)]; 
f1_hat = mu; 
f2_hat = zeros(4,1);
tau_prev = zeros(4,1); 
rho_hat_filt = zeros(4,1);

log_eta = zeros(4, N); log_rho = zeros(4, N); log_rho_hat = zeros(4, N);
log_t = zeros(1, N); log_tau = zeros(3, N);
log_psi = zeros(1, N); log_psi_d = zeros(1, N);

R_accept = 1.5; 
tau_rho_filt = 0.2;
D_nom_val = diag([params.d11, params.d22, params.d33, params.d44]);
wd_filt = 0; 
tau_w_filt = 1.0;

% ============================================================
% 主循环
% ============================================================
for k = 1:N
    tk = (k-1)*params.Ts;
    log_eta(:,k) = eta; 
    log_t(k) = tk;
    
    % ---------- LOS 制导 ----------
    [psi_d, psi_continuous, target_idx] = ...
        los_guidance(eta, mu, path, target_idx, psi_continuous, R_accept, Delta);
    log_psi(k) = eta(4); 
    log_psi_d(k) = psi_d;
    
    % ---------- 垂向制导 ----------
    dists_xy = sqrt((path(:,1) - eta(1)).^2 + (path(:,2) - eta(2)).^2);
    [~, nearest_idx] = min(dists_xy);
    z_ref = path(nearest_idx, 3);
    wd_raw = params.kz * (z_ref - eta(3));
    wd_raw = max(min(wd_raw, params.wmax), -params.wmax);
    wd_filt = wd_filt + params.Ts/tau_w_filt * (wd_raw - wd_filt);
    wd = wd_filt;

    if tk < 60
    ud = 0.05;    % 前 30 s 降速到 0.06 m/s
else
    ud = 0.1;    % 之后恢复正常
end
    
    % ---------- 航向误差 → r_d ----------
    e_psi = psi_d - eta(4);
    while e_psi > pi,  e_psi = e_psi - 2*pi; end
    while e_psi < -pi, e_psi = e_psi + 2*pi; end
    r_d = 0.15 * e_psi;
    r_d = max(min(r_d, 0.1), -0.1);
    mu_d = [ud; 0; wd; r_d];
    
    % ---------- ESO 扰动估计 ----------
    [~, rho_hat, f1_hat, f2_hat] = ...
        eso_estimator(mu, tau_prev, f1_hat, f2_hat, params);
    log_rho_hat(:,k) = rho_hat;
    
    % ---------- MPC ----------
    [tau_bar_nom, xi] = mpc_controller(mu, mu_d, xi, params);
    
    % ---------- 前馈补偿（选择性：仅 u、v）----------
    tau_bar_full = [tau_bar_nom(1); 0; tau_bar_nom(2); tau_bar_nom(3)];
    rho_hat_filt = rho_hat_filt + params.Ts/tau_rho_filt * (rho_hat - rho_hat_filt);
    rho_nominal = -params.Minv_nom * D_nom_val * mu;
    rho_hat_extra = rho_hat_filt - rho_nominal;
    soft_start = min(tk / 3.0, 1.0);
    rho_hat_comp = soft_start * rho_hat_extra;
    rho_hat_comp = max(min(rho_hat_comp, 0.05), -0.05);
    rho_hat_comp(3) = 0;
    rho_hat_comp(4) = 0;
    tau_full = tau_bar_full - params.M_nom * rho_hat_comp;
    tau = [tau_full(1); 0; tau_full(3); tau_full(4)];
    
    % ---------- 限幅 ----------
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
    
    % 提前到达终点则退出
    if norm(eta(1:3) - path(end,:)') < 1.0
        fprintf('AUV 到达终点，t = %.2f s\n', tk);
        break;
    end
end

% 截断到实际步数
N_used = k;
log_eta = log_eta(:, 1:N_used);
log_rho = log_rho(:, 1:N_used);
log_rho_hat = log_rho_hat(:, 1:N_used);
log_t = log_t(1:N_used);
log_tau = log_tau(:, 1:N_used);
log_psi = log_psi(1:N_used);
log_psi_d = log_psi_d(1:N_used);

% ============================================================
% 计算 RMS
% ============================================================
errors_x = zeros(1, N_used); 
errors_y = zeros(1, N_used); 
errors_z = zeros(1, N_used);
for k = 1:N_used
    [~, idx] = min(abs(path(:,1) - log_eta(1,k)));
    errors_x(k) = abs(path(idx,1) - log_eta(1,k));
    errors_y(k) = abs(path(idx,2) - log_eta(2,k));
    errors_z(k) = abs(path(idx,3) - log_eta(3,k));
end
errors_3d = sqrt(errors_x.^2 + errors_y.^2 + errors_z.^2);
rms_x = sqrt(mean(errors_x.^2));
rms_y = sqrt(mean(errors_y.^2));
rms_z = sqrt(mean(errors_z.^2));
rms_3d = sqrt(mean(errors_3d.^2));
max_err = max(errors_3d);

fprintf('==========================================\n');
fprintf('A*+APF 路径跟踪结果\n');
fprintf('X 方向 RMS: %.4f m\n', rms_x);
fprintf('Y 方向 RMS: %.4f m\n', rms_y);
fprintf('Z 方向 RMS: %.4f m\n', rms_z);
fprintf('三维合成 RMS: %.4f m\n', rms_3d);
fprintf('最大误差: %.4f m\n', max_err);
fprintf('==========================================\n');


% ================================================================
% 图 4.7 三维轨迹（含障碍物、起点、终点）
% ================================================================
% 加载障碍物数据（与规划代码保持一致）
load('astar_path.mat', 'path', 'start_pos', 'goal_pos');
% 手动定义障碍物（与规划代码中的obs一致）
obstacles = [
    5, 5, 5;
    9, 7.5, 9.5;
    7, 8, 9;
    8, 8, 8;
    10, 10, 10;
    15, 15, 15;
    18, 18, 20;
    20, 10, 25;
    25, 25, 25;
    15, 15, 15   % 中心点
];

figure('Position', [100, 100, 900, 700]);
hold on;

% 绘制静态障碍物
for o = 1:size(obstacles,1)
    scatter3(obstacles(o,1), obstacles(o,2), obstacles(o,3), ...
             200, 'k', 'filled', 'DisplayName', '障碍物');
end

% 绘制规划路径
h1 = plot3(path(:,1), path(:,2), path(:,3), 'k--', 'LineWidth', 1.8, ...
           'DisplayName', 'A*+APF 规划路径');

% 绘制跟踪轨迹
h2 = plot3(log_eta(1,:), log_eta(2,:), log_eta(3,:), 'b-', 'LineWidth', 1.5, ...
           'DisplayName', 'ESO-MPC 跟踪');

% 绘制起点和终点
h3 = plot3(path(1,1), path(1,2), path(1,3), 'g*', ...
           'MarkerSize', 18, 'LineWidth', 2, 'DisplayName', '起点');
h4 = plot3(path(end,1), path(end,2), path(end,3), 'r*', ...
           'MarkerSize', 18, 'LineWidth', 2, 'DisplayName', '终点');

xlabel('X (m)', 'FontSize', 12);
ylabel('Y (m)', 'FontSize', 12);
zlabel('Z (m)', 'FontSize', 12);
legend([h1, h2, h3, h4], 'FontSize', 11, 'Location', 'best');
grid on; view(45, 30);
title('A*+APF 规划路径跟踪', 'FontSize', 14);
axis([0 30 0 30 0 30]);
print('Fig4_7', '-dpng', '-r300');


% ================================================================
% 图 4.8 XY 和 XZ 投影（含障碍物）
% ================================================================
figure('Position', [100, 100, 1000, 500]);

% --- XY 投影 ---
subplot(1,2,1);
hold on;
% 障碍物投影（用矩形框表示）
for o = 1:size(obstacles,1)
    plot(obstacles(o,1), obstacles(o,2), 'ks', ...
         'MarkerSize', 12, 'MarkerFaceColor', 'k');
end
h1 = plot(path(:,1), path(:,2), 'k--', 'LineWidth', 1.5, 'DisplayName', '规划');
h2 = plot(log_eta(1,:), log_eta(2,:), 'b-', 'LineWidth', 1.2, 'DisplayName', 'ESO-MPC');
h3 = plot(path(1,1), path(1,2), 'g*', 'MarkerSize', 15, 'LineWidth', 2, ...
          'DisplayName', '起点');
h4 = plot(path(end,1), path(end,2), 'r*', 'MarkerSize', 15, 'LineWidth', 2, ...
          'DisplayName', '终点');
xlabel('X (m)', 'FontSize', 12); ylabel('Y (m)', 'FontSize', 12);
legend([h1,h2,h3,h4], 'FontSize', 10, 'Location', 'best');
grid on; axis equal;
title('XY 投影', 'FontSize', 13);
xlim([0 30]); ylim([0 30]);

% --- XZ 投影 ---
subplot(1,2,2);
hold on;
for o = 1:size(obstacles,1)
    plot(obstacles(o,1), obstacles(o,3), 'ks', ...
         'MarkerSize', 12, 'MarkerFaceColor', 'k');
end
plot(path(:,1), path(:,3), 'k--', 'LineWidth', 1.5);
plot(log_eta(1,:), log_eta(3,:), 'b-', 'LineWidth', 1.2);
plot(path(1,1), path(1,3), 'g*', 'MarkerSize', 15, 'LineWidth', 2);
plot(path(end,1), path(end,3), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
xlabel('X (m)', 'FontSize', 12); ylabel('Z (m)', 'FontSize', 12);
grid on; title('XZ 投影', 'FontSize', 13);
xlim([0 30]); ylim([0 30]);

print('Fig4_8', '-dpng', '-r300');

% ================================================================
% 图 4.9 三方向误差
% ================================================================
figure('Position', [100, 100, 800, 600]);
subplot(3,1,1);
plot(log_t, errors_x, 'b-', 'LineWidth', 1.2);
ylabel('E_x (m)', 'FontSize', 12); grid on;
title('A*+APF 规划路径跟踪误差', 'FontSize', 13);
subplot(3,1,2);
plot(log_t, errors_y, 'r-', 'LineWidth', 1.2);
ylabel('E_y (m)', 'FontSize', 12); grid on;
subplot(3,1,3);
plot(log_t, errors_z, 'g-', 'LineWidth', 1.2);
xlabel('时间 (s)', 'FontSize', 12); 
ylabel('E_z (m)', 'FontSize', 12);
grid on;
print('Fig4_9', '-dpng', '-r300');

% ================================================================
% 图 4.10 控制输入
% ================================================================
figure('Position', [100, 100, 800, 500]);
plot(log_t, log_tau(1,:), 'b-', 'LineWidth', 1.2); hold on;
plot(log_t, log_tau(2,:), 'r-', 'LineWidth', 1.2);
plot(log_t, log_tau(3,:), 'g-', 'LineWidth', 1.2);
xlabel('时间 (s)', 'FontSize', 12);
ylabel('控制量 (N, N·m)', 'FontSize', 12);
legend('\tau_u', '\tau_w', '\tau_r', 'FontSize', 11);
grid on; title('A*+APF 路径控制输入', 'FontSize', 13);
print('Fig4_10', '-dpng', '-r300');

% ================================================================
% 图 4.11 ESO 扰动估计
% ================================================================
channels = {'u', 'v', 'w', 'r'};
figure('Position', [100, 100, 800, 700]);
for i = 1:4
    subplot(4,1,i);
    plot(log_t, log_rho(i,:), 'k--', 'LineWidth', 1.2); hold on;
    plot(log_t, log_rho_hat(i,:), 'r-', 'LineWidth', 1.2);
    ylabel(['\rho_', channels{i}], 'FontSize', 11); grid on;
    if i == 1
        legend('真实扰动', 'ESO估计', 'FontSize', 10, 'Location', 'best');
        title('ESO 扰动估计效果', 'FontSize', 13);
    end
    if i == 4, xlabel('时间 (s)', 'FontSize', 11); end
end
print('Fig4_11', '-dpng', '-r300');

save('result_4_3.mat', 'path', 'log_eta', 'log_t', 'log_tau', ...
     'log_rho', 'log_rho_hat', 'log_psi', 'log_psi_d', ...
     'rms_x', 'rms_y', 'rms_z', 'rms_3d', 'max_err');