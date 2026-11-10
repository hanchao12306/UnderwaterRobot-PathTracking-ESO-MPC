clear; clc; close all;
clear functions;
utils_plot();
params = set_params();

% ============================================================
% 加载 A*+APF 规划路径（与 test9 一致）
% ============================================================
if exist('astar_path.mat', 'file')
    load('astar_path.mat', 'path', 'start_pos', 'goal_pos');
    fprintf('已加载 A*+APF 规划路径，点数: %d\n', size(path,1));
    fprintf('起点: (%.2f, %.2f, %.2f)\n', start_pos(1), start_pos(2), start_pos(3));
    fprintf('终点: (%.2f, %.2f, %.2f)\n', goal_pos(1), goal_pos(2), goal_pos(3));
else
    error('未找到 astar_path.mat，请先运行规划代码并保存');
end

if size(path, 1) < 200
    t = 1:size(path,1);
    tq = linspace(1, size(path,1), 300);
    path = [interp1(t, path(:,1), tq, 'spline')', ...
            interp1(t, path(:,2), tq, 'spline')', ...
            interp1(t, path(:,3), tq, 'spline')'];
end

% 计算路径总长度，估算仿真时长（与 test9 一致）
path_len = 0;
for i = 1:size(path,1)-1
    path_len = path_len + norm(path(i+1,:) - path(i,:));
end
fprintf('路径总长度: %.2f m\n', path_len);

ud_ref = 0.1;
T_est = path_len / ud_ref + 20;
N = ceil(T_est / params.Ts);
fprintf('预计仿真时长: %.1f s，步数: %d\n', T_est, N);

% ============================================================
% 障碍物列表（与规划代码完全一致，共 10 个）
% ============================================================
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
    15, 15, 15
];

% 起点/终点（整数，与规划图一致）
start_node_plot = [3, 4, 5];
target_node_plot = [26, 27, 28];

% 动态障碍物初始位置
moving_obs_init = [9, 11, 12];

% ============================================================
% 运行三种方案
% ============================================================
fprintf('\n正在运行 ESO-MPC ...\n');
[eta1, tau1] = run_esompc(path, params, N);

fprintf('正在运行 标准MPC ...\n');
[eta2, tau2] = run_mpc(path, params, N);

fprintf('正在运行 PID ...\n');
[eta3, tau3] = run_pid(path, params, N);

% ============================================================
% 计算 RMS
% ============================================================
[rms1, max1, ex1, ey1, ez1, t1] = compute_rms_3d(eta1, path, params);
[rms2, max2, ex2, ey2, ez2, t2] = compute_rms_3d(eta2, path, params);
[rms3, max3, ex3, ey3, ez3, t3] = compute_rms_3d(eta3, path, params);

fprintf('\n==========================================\n');
fprintf('      方法        三维RMS   最大误差\n');
fprintf('    ESO-MPC      %.4f    %.4f\n', rms1, max1);
fprintf('    标准 MPC     %.4f    %.4f\n', rms2, max2);
fprintf('    PID          %.4f    %.4f\n', rms3, max3);
fprintf('==========================================\n');
fprintf('ESO-MPC 相比标准 MPC 提升: %.1f%%\n', (rms2-rms1)/rms2*100);
fprintf('ESO-MPC 相比 PID 提升:     %.1f%%\n', (rms3-rms1)/rms3*100);
fprintf('==========================================\n');

% ============================================================
% 图 5.7 三维轨迹对比（图例只显示 4 条路径）
% ============================================================
figure('Name', 'A-APF 三维路径规划结果', 'Position', [80, 100, 600, 600]);
view(3); grid on; hold on;
axis equal;
xlim([0, 30]); ylim([0, 30]); zlim([0, 30]);
xlabel('X轴', 'FontSize', 11);
ylabel('Y轴', 'FontSize', 11);
zlabel('Z轴', 'FontSize', 11);

% ---------- 起点（不进图例） ----------
plot3(start_node_plot(1), start_node_plot(2), start_node_plot(3), ...
      'g*', 'MarkerSize', 10);

% ---------- 终点（不进图例） ----------
plot3(target_node_plot(1), target_node_plot(2), target_node_plot(3), ...
      'r*', 'MarkerSize', 10);

% ---------- 静态障碍物（不进图例） ----------
scatter3(obstacles(:,1), obstacles(:,2), obstacles(:,3), ...
         150, 'k', 'filled');

% ---------- 动态障碍物球体（不进图例） ----------
[x_sphere, y_sphere, z_sphere] = sphere(20);
radius = 0.5;
x_sphere = x_sphere * radius + moving_obs_init(1);
y_sphere = y_sphere * radius + moving_obs_init(2);
z_sphere = z_sphere * radius + moving_obs_init(3);
surf(x_sphere, y_sphere, z_sphere, ...
     'FaceColor', 'cyan', 'FaceAlpha', 0.8, 'EdgeColor', 'none', ...
     'HandleVisibility', 'off');

% ---------- 4 条路径（进图例） ----------
h_plan = plot3(path(:,1), path(:,2), path(:,3), 'k--', 'LineWidth', 1.8, ...
               'DisplayName', 'A*+APF 规划路径');
h_eso  = plot3(eta1(1,:), eta1(2,:), eta1(3,:), 'b-', 'LineWidth', 1.5, ...
               'DisplayName', 'ESO-MPC');
h_mpc  = plot3(eta2(1,:), eta2(2,:), eta2(3,:), 'r-', 'LineWidth', 1.5, ...
               'DisplayName', '标准MPC');
h_pid  = plot3(eta3(1,:), eta3(2,:), eta3(3,:), 'g-', 'LineWidth', 1.5, ...
               'DisplayName', 'PID');

% ---------- 图例：只列 4 条路径 ----------
legend([h_plan, h_eso, h_mpc, h_pid], ...
       {'A*+APF 规划路径', 'ESO-MPC', '标准MPC', 'PID'}, ...
       'Location', 'best', 'FontSize', 10);

title('A*+APF 规划路径跟踪对比', 'FontSize', 13);
print('Fig5_7', '-dpng', '-r300');
% ============================================================
% 【新增】图4.13放大子图：障碍物密集区
% ============================================================
ax_main = gca;
pos_main = get(ax_main, 'Position');

% ---------- 创建内嵌放大子图（右上角）----------
ax_zoom = axes('Position', ...
    [pos_main(1) + pos_main(3)*0.52, ...    % 水平位置
     pos_main(2) + pos_main(4)*0.65, ...    % 垂直位置
     pos_main(3)*0.40, ...                  % 宽度
     pos_main(4)*0.40]);                    % 高度
hold(ax_zoom, 'on');
box(ax_zoom, 'on');

% ---------- 在放大子图中重绘数据 ----------
plot3(ax_zoom, path(:,1), path(:,2), path(:,3), 'k--', 'LineWidth', 1.5);
plot3(ax_zoom, eta1(1,:), eta1(2,:), eta1(3,:), 'b-', 'LineWidth', 1.2);
plot3(ax_zoom, eta2(1,:), eta2(2,:), eta2(3,:), 'r-', 'LineWidth', 1.2);
plot3(ax_zoom, eta3(1,:), eta3(2,:), eta3(3,:), 'g-', 'LineWidth', 1.2);

% 绘制该区域障碍物
for o = 1:size(obstacles,1)
    scatter3(ax_zoom, obstacles(o,1), obstacles(o,2), obstacles(o,3), ...
             80, 'k', 'filled');
end

% 起点终点
plot3(ax_zoom, start_node_plot(1), start_node_plot(2), start_node_plot(3), ...
      'g*', 'MarkerSize', 12);
plot3(ax_zoom, target_node_plot(1), target_node_plot(2), target_node_plot(3), ...
      'r*', 'MarkerSize', 12);

% ---------- 设置放大区域（障碍物密集区 3~12 m）----------
xlim(ax_zoom, [3 12]);
ylim(ax_zoom, [3 12]);
zlim(ax_zoom, [3 12]);
view(ax_zoom, 45, 30);
grid(ax_zoom, 'on');
set(ax_zoom, 'FontSize', 8);
title(ax_zoom, '障碍物密集区放大', 'FontSize', 9);

% 主图上绘制放大范围指示框（灰色虚线，8条棱）
% 底面
plot3(ax_main, [3 12 12 3 3], [3 3 12 12 3], [3 3 3 3 3], ...
      'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'LineStyle', '--', ...
      'HandleVisibility', 'off');
% 顶面
plot3(ax_main, [3 12 12 3 3], [3 3 12 12 3], [12 12 12 12 12], ...
      'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'LineStyle', '--', ...
      'HandleVisibility', 'off');
% 4条竖棱
for xy = [[3 3]; [12 3]; [12 12]; [3 12]]'
    plot3(ax_main, [xy(1) xy(1)], [xy(2) xy(2)], [3 12], ...
          'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'LineStyle', '--', ...
          'HandleVisibility', 'off');
end
% ============================================================

% ============================================================
% 图 5.8 XY/XZ 投影对比（两子图比例一致）
% ============================================================
figure('Position', [100, 100, 1000, 500]);

% --- XY 投影 ---
subplot(1,2,1);
hold on;
plot(obstacles(:,1), obstacles(:,2), 'ks', ...
     'MarkerSize', 12, 'MarkerFaceColor', 'k');
plot(path(:,1), path(:,2), 'k--', 'LineWidth', 1.5);
plot(eta1(1,:), eta1(2,:), 'b-', 'LineWidth', 1.2);
plot(eta2(1,:), eta2(2,:), 'r-', 'LineWidth', 1.2);
plot(eta3(1,:), eta3(2,:), 'g-', 'LineWidth', 1.2);
plot(start_node_plot(1), start_node_plot(2), 'g*', ...
     'MarkerSize', 12, 'LineWidth', 2);
plot(target_node_plot(1), target_node_plot(2), 'r*', ...
     'MarkerSize', 12, 'LineWidth', 2);
xlabel('X (m)', 'FontSize', 12);
ylabel('Y (m)', 'FontSize', 12);
legend('障碍物', '规划', 'ESO-MPC', 'MPC', 'PID', '起点', '终点', ...
       'FontSize', 9, 'Location', 'best');
grid on;
axis equal;                      % ✅ XY 等比例
xlim([0 30]); ylim([0 30]);
title('XY 投影', 'FontSize', 13);

% --- XZ 投影 ---
subplot(1,2,2);
hold on;
plot(obstacles(:,1), obstacles(:,3), 'ks', ...
     'MarkerSize', 12, 'MarkerFaceColor', 'k');
plot(path(:,1), path(:,3), 'k--', 'LineWidth', 1.5);
plot(eta1(1,:), eta1(3,:), 'b-', 'LineWidth', 1.2);
plot(eta2(1,:), eta2(3,:), 'r-', 'LineWidth', 1.2);
plot(eta3(1,:), eta3(3,:), 'g-', 'LineWidth', 1.2);
plot(start_node_plot(1), start_node_plot(3), 'g*', ...
     'MarkerSize', 12, 'LineWidth', 2);
plot(target_node_plot(1), target_node_plot(3), 'r*', ...
     'MarkerSize', 12, 'LineWidth', 2);
xlabel('X (m)', 'FontSize', 12);
ylabel('Z (m)', 'FontSize', 12);
grid on;
axis equal;                      % ✅ XZ 也加等比例
xlim([0 30]); ylim([0 30]);
title('XZ 投影', 'FontSize', 13);

print('Fig5_8', '-dpng', '-r300');

% ============================================================
% 图 5.9 三方向误差对比
% ============================================================
figure('Position', [100, 100, 900, 800]);
subplot(3,1,1);
plot(t1, ex1, 'b-', 'LineWidth', 1.2); hold on;
plot(t2, ex2, 'r--', 'LineWidth', 1.2);
plot(t3, ex3, 'g-.', 'LineWidth', 1.2);
ylabel('E_x (m)', 'FontSize', 12); grid on;
legend('ESO-MPC', '标准MPC', 'PID', 'FontSize', 10);
title('三维路径跟踪误差对比', 'FontSize', 13);

subplot(3,1,2);
plot(t1, ey1, 'b-', 'LineWidth', 1.2); hold on;
plot(t2, ey2, 'r--', 'LineWidth', 1.2);
plot(t3, ey3, 'g-.', 'LineWidth', 1.2);
ylabel('E_y (m)', 'FontSize', 12); grid on;

subplot(3,1,3);
plot(t1, ez1, 'b-', 'LineWidth', 1.2); hold on;
plot(t2, ez2, 'r--', 'LineWidth', 1.2);
plot(t3, ez3, 'g-.', 'LineWidth', 1.2);
xlabel('时间 (s)', 'FontSize', 12);
ylabel('E_z (m)', 'FontSize', 12); grid on;
print('Fig5_9', '-dpng', '-r300');
% ============================================================
% 【新增】图4.18放大子图：Z方向误差放大
% ============================================================
ax_main = gca;
pos_main = get(ax_main, 'Position');

% ---------- 创建内嵌放大子图（右上角）----------
ax_zoom = axes('Position', ...
    [pos_main(1) + pos_main(3)*0.52, ...    % 水平位置
     pos_main(2) + pos_main(4)*0.45, ...    % 垂直位置
     pos_main(3)*0.36, ...                  % 宽度
     pos_main(4)*0.45]);                    % 高度
hold(ax_zoom, 'on');
box(ax_zoom, 'on');

% ---------- 重绘Z方向误差 ----------
plot(ax_zoom, t1, ez1, 'b-', 'LineWidth', 1.2);
plot(ax_zoom, t2, ez2, 'r--', 'LineWidth', 1.2);
plot(ax_zoom, t3, ez3, 'g-.', 'LineWidth', 1.2);

% ---------- 设置放大区域（t=150~350 s，PID误差累积段）----------
xlim(ax_zoom, [150 350]);
ylim(ax_zoom, [0 5]);
grid(ax_zoom, 'on');
set(ax_zoom, 'FontSize', 8);
title(ax_zoom, 'Z方向误差放大', 'FontSize', 9);
xlabel(ax_zoom, 't (s)', 'FontSize', 8);

% ---------- 在主图上标注放大范围（红色虚线框）----------
rectangle('Position', [150, 0, 200, 5], ...
          'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
% ============================================================

% ============================================================
% 图 5.10 RMS误差柱状图
% ============================================================
figure('Position', [100, 100, 700, 500]);
rms_data = [rms1, rms2, rms3];
bar(rms_data, 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTickLabel', {'ESO-MPC', '标准MPC', 'PID'});
ylabel('三维 RMS 误差 (m)', 'FontSize', 12);
grid on;
title('三种方案 RMS 误差对比', 'FontSize', 13);
for i = 1:3
    text(i, rms_data(i)+0.02, sprintf('%.3f', rms_data(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
end
print('Fig5_10', '-dpng', '-r300');

save('result_compare.mat', 'path', 'eta1', 'eta2', 'eta3', ...
     'rms1', 'rms2', 'rms3', 'max1', 'max2', 'max3');
fprintf('\n所有图表已保存 (Fig5_7 ~ Fig5_10)\n');


% ============================================================
% 子函数1：ESO-MPC
% ============================================================
function [eta_log, tau_log] = run_esompc(path, params, N)
    eta = [path(1,1); path(1,2); path(1,3); 0];
    psi0 = atan2(path(5,2)-path(1,2), path(5,1)-path(1,1));
    eta(4) = psi0;
    mu = [0.3; 0; 0; 0];
    psi_continuous = psi0; target_idx = 1; Delta = 0.8;
    xi = [mu; zeros(3,1)]; f1_hat = mu; f2_hat = zeros(4,1);
    tau_prev = zeros(4,1); rho_hat_filt = zeros(4,1);
    eta_log = zeros(4, N); tau_log = zeros(3, N);
    R_accept = 1.5; tau_rho_filt = 0.2;
    D_nom_val = diag([params.d11, params.d22, params.d33, params.d44]);
    wd_filt = 0; tau_w_filt = 1.0;
    
    for k = 1:N
        tk = (k-1)*params.Ts;
        eta_log(:,k) = eta;
        [psi_d, psi_continuous, target_idx] = ...
            los_guidance(eta, mu, path, target_idx, psi_continuous, R_accept, Delta);
        dists_xy = sqrt((path(:,1)-eta(1)).^2 + (path(:,2)-eta(2)).^2);
        [~, nearest_idx] = min(dists_xy);
        z_ref = path(nearest_idx, 3);
        wd_raw = params.kz * (z_ref - eta(3));
        wd_raw = max(min(wd_raw, params.wmax), -params.wmax);
        wd_filt = wd_filt + params.Ts/tau_w_filt * (wd_raw - wd_filt);
        wd = wd_filt;
        if tk < 60, ud = 0.05; else, ud = 0.1; end
        e_psi = psi_d - eta(4);
        while e_psi > pi,  e_psi = e_psi - 2*pi; end
        while e_psi < -pi, e_psi = e_psi + 2*pi; end
        r_d = 0.15 * e_psi;
        r_d = max(min(r_d, 0.1), -0.1);
        mu_d = [ud; 0; wd; r_d];
        [~, rho_hat, f1_hat, f2_hat] = ...
            eso_estimator(mu, tau_prev, f1_hat, f2_hat, params);
        [tau_bar_nom, xi] = mpc_controller(mu, mu_d, xi, params);
        tau_bar_full = [tau_bar_nom(1); 0; tau_bar_nom(2); tau_bar_nom(3)];
        rho_hat_filt = rho_hat_filt + params.Ts/tau_rho_filt * (rho_hat - rho_hat_filt);
        rho_nominal = -params.Minv_nom * D_nom_val * mu;
        rho_hat_extra = rho_hat_filt - rho_nominal;
        soft_start = min(tk/3.0, 1.0);
        rho_hat_comp = soft_start * rho_hat_extra;
        rho_hat_comp = max(min(rho_hat_comp, 0.05), -0.05);
        rho_hat_comp(3) = 0; rho_hat_comp(4) = 0;
        tau_full = tau_bar_full - params.M_nom * rho_hat_comp;
        tau = [tau_full(1); 0; tau_full(3); tau_full(4)];
        tau(1) = max(min(tau(1), params.tau_max(1)), params.tau_min(1));
        tau(3) = max(min(tau(3), params.tau_max(2)), params.tau_min(2));
        tau(4) = max(min(tau(4), params.tau_max(3)), params.tau_min(3));
        tau_prev = tau; tau_log(:,k) = tau([1,3,4]);
        [mu, ~] = auv_dynamics(mu, tau, tk, params);
        psi = eta(4);
        J = [cos(psi), -sin(psi), 0, 0; sin(psi), cos(psi), 0, 0;
             0, 0, 1, 0; 0, 0, 0, 1];
        eta = eta + params.Ts * J * mu;
        if norm(eta(1:3) - path(end,:)') < 1.0
            eta_log = eta_log(:, 1:k); tau_log = tau_log(:, 1:k); return;
        end
    end
end


% ============================================================
% 子函数2：标准 MPC
% ============================================================
function [eta_log, tau_log] = run_mpc(path, params, N)
    eta = [path(1,1); path(1,2); path(1,3); 0];
    psi0 = atan2(path(5,2)-path(1,2), path(5,1)-path(1,1));
    eta(4) = psi0; mu = [0.3; 0; 0; 0];
    psi_continuous = psi0; target_idx = 1; Delta = 0.8;
    xi = [mu; zeros(3,1)];
    eta_log = zeros(4, N); tau_log = zeros(3, N);
    R_accept = 1.5; wd_filt = 0; tau_w_filt = 1.0;
    for k = 1:N
        tk = (k-1)*params.Ts;
        eta_log(:,k) = eta;
        [psi_d, psi_continuous, target_idx] = ...
            los_guidance(eta, mu, path, target_idx, psi_continuous, R_accept, Delta);
        dists_xy = sqrt((path(:,1)-eta(1)).^2 + (path(:,2)-eta(2)).^2);
        [~, nearest_idx] = min(dists_xy);
        z_ref = path(nearest_idx, 3);
        wd_raw = params.kz * (z_ref - eta(3));
        wd_raw = max(min(wd_raw, params.wmax), -params.wmax);
        wd_filt = wd_filt + params.Ts/tau_w_filt * (wd_raw - wd_filt);
        wd = wd_filt;
        if tk < 60, ud = 0.05; else, ud = 0.1; end
        e_psi = psi_d - eta(4);
        while e_psi > pi, e_psi = e_psi - 2*pi; end
        while e_psi < -pi, e_psi = e_psi + 2*pi; end
        r_d = 0.15 * e_psi; r_d = max(min(r_d, 0.1), -0.1);
        mu_d = [ud; 0; wd; r_d];
        [tau_bar_nom, xi] = mpc_controller(mu, mu_d, xi, params);
        tau_bar_full = [tau_bar_nom(1); 0; tau_bar_nom(2); tau_bar_nom(3)];
        tau = [tau_bar_full(1); 0; tau_bar_full(3); tau_bar_full(4)];
        tau(1) = max(min(tau(1), params.tau_max(1)), params.tau_min(1));
        tau(3) = max(min(tau(3), params.tau_max(2)), params.tau_min(2));
        tau(4) = max(min(tau(4), params.tau_max(3)), params.tau_min(3));
        tau_log(:,k) = tau([1,3,4]);
        [mu, ~] = auv_dynamics(mu, tau, tk, params);
        psi = eta(4);
        J = [cos(psi), -sin(psi), 0, 0; sin(psi), cos(psi), 0, 0;
             0, 0, 1, 0; 0, 0, 0, 1];
        eta = eta + params.Ts * J * mu;
        if norm(eta(1:3) - path(end,:)') < 1.0
            eta_log = eta_log(:, 1:k); tau_log = tau_log(:, 1:k); return;
        end
    end
end


% ============================================================
% 子函数3：PID
% ============================================================
function [eta_log, tau_log] = run_pid(path, params, N)
    eta = [path(1,1); path(1,2); path(1,3); 0];
    psi0 = atan2(path(5,2)-path(1,2), path(5,1)-path(1,1));
    eta(4) = psi0; mu = [0.3; 0; 0; 0];
    psi_continuous = psi0; target_idx = 1; Delta = 0.8;
    eta_log = zeros(4, N); tau_log = zeros(3, N);
    R_accept = 1.5; wd_filt = 0; tau_w_filt = 1.0;
    integral_psi = 0; integral_u = 0;
    for k = 1:N
        tk = (k-1)*params.Ts;
        eta_log(:,k) = eta;
        [psi_d, psi_continuous, target_idx] = ...
            los_guidance(eta, mu, path, target_idx, psi_continuous, R_accept, Delta);
        dists_xy = sqrt((path(:,1)-eta(1)).^2 + (path(:,2)-eta(2)).^2);
        [~, nearest_idx] = min(dists_xy);
        z_ref = path(nearest_idx, 3);
        wd_raw = params.kz * (z_ref - eta(3));
        wd_raw = max(min(wd_raw, params.wmax), -params.wmax);
        wd_filt = wd_filt + params.Ts/tau_w_filt * (wd_raw - wd_filt);
        wd = wd_filt;
        if tk < 60, ud = 0.05; else, ud = 0.1; end
        e_psi = psi_d - eta(4);
        while e_psi > pi, e_psi = e_psi - 2*pi; end
        while e_psi < -pi, e_psi = e_psi + 2*pi; end
        e_u = ud - mu(1);
        Kp_psi = 20; Ki_psi = 0.5; Kd_psi = 15;
        Kp_u = 40; Ki_u = 1;
        integral_psi = max(min(integral_psi + e_psi*params.Ts, 1), -1);
        integral_u = max(min(integral_u + e_u*params.Ts, 1), -1);
        tau_r = Kp_psi*e_psi + Ki_psi*integral_psi - Kd_psi*mu(4);
        tau_u = Kp_u*e_u + Ki_u*integral_u;
        tau_w = 30 * (wd - mu(3));
        tau_u = max(min(tau_u, params.tau_max(1)), params.tau_min(1));
        tau_w = max(min(tau_w, params.tau_max(2)), params.tau_min(2));
        tau_r = max(min(tau_r, params.tau_max(3)), params.tau_min(3));
        tau = [tau_u; 0; tau_w; tau_r];
        tau_log(:,k) = tau([1,3,4]);
        [mu, ~] = auv_dynamics(mu, tau, tk, params);
        psi = eta(4);
        J = [cos(psi), -sin(psi), 0, 0; sin(psi), cos(psi), 0, 0;
             0, 0, 1, 0; 0, 0, 0, 1];
        eta = eta + params.Ts * J * mu;
        if norm(eta(1:3) - path(end,:)') < 1.0
            eta_log = eta_log(:, 1:k); tau_log = tau_log(:, 1:k); return;
        end
    end
end


% ============================================================
% 子函数4：RMS计算
% ============================================================
function [rms_3d, max_err, ex, ey, ez, t_log] = compute_rms_3d(eta_log, path, params)
    N_used = size(eta_log, 2);
    ex = zeros(1, N_used); ey = zeros(1, N_used); ez = zeros(1, N_used);
    for k = 1:N_used
        [~, idx] = min(abs(path(:,1) - eta_log(1,k)));
        ex(k) = abs(path(idx,1) - eta_log(1,k));
        ey(k) = abs(path(idx,2) - eta_log(2,k));
        ez(k) = abs(path(idx,3) - eta_log(3,k));
    end
    err_3d = sqrt(ex.^2 + ey.^2 + ez.^2);
    rms_3d = sqrt(mean(err_3d.^2));
    max_err = max(err_3d);
    t_log = (0:N_used-1) * params.Ts;
end

% ============================================================
% 计算综合评价指标（最终定稿版）
% ============================================================
fprintf('\n=== 综合评价指标 ===\n');
fprintf('%-10s %-10s %-10s %-10s %-10s\n', ...
        '方法', 'RMS(m)', '最大(m)', '稳态(m)', '收敛性');
fprintf('%s\n', repmat('-', 1, 55));

for m = 1:3
    if m == 1
        eta = eta1; name = 'ESO-MPC';
    elseif m == 2
        eta = eta2; name = '标准MPC';
    else
        eta = eta3; name = 'PID';
    end
    
    N_used = size(eta, 2);
    err_3d = zeros(1, N_used);
    for k = 1:N_used
        [~, idx] = min(abs(path(:,1) - eta(1,k)));
        err_3d(k) = sqrt((path(idx,1)-eta(1,k))^2 + ...
                         (path(idx,2)-eta(2,k))^2 + ...
                         (path(idx,3)-eta(3,k))^2);
    end
    
    rms_val = sqrt(mean(err_3d.^2));
    max_val = max(err_3d);
    
    % 稳态误差（最后30%时间）
    steady_start = round(N_used * 0.7);
    steady_val = sqrt(mean(err_3d(steady_start:end).^2));
    
    % 收敛性判定：稳态误差 < 1.0 m 视为收敛
    if steady_val < 1.0
        conv_str = '收敛';
    else
        conv_str = '未收敛';
    end
    
    fprintf('%-10s %-10.4f %-10.4f %-10.4f %-10s\n', ...
            name, rms_val, max_val, steady_val, conv_str);
end
fprintf('\n');