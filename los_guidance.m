function [psi_d, psi_continuous, target_idx] = ...
    los_guidance(eta, mu, path, target_idx, psi_continuous, R_accept, Delta)
% ============================================================
% LOS 制导（单调递增最近点搜索）
% 用传入的 target_idx 作为搜索起点，只允许向前搜索
% 避免圆路径起点终点重合导致的跳变
% ============================================================

x = eta(1); y = eta(2); z = eta(3);
u = mu(1); v = mu(2);

% 【核心】从上一次 target_idx 开始，向前搜索 40 个点
search_start = max(1, target_idx - 5);         % 允许小幅回退
search_end   = min(target_idx + 40, size(path,1));
search_range = search_start:search_end;

dists = sqrt((path(search_range,1) - x).^2 + ...
             (path(search_range,2) - y).^2);
[~, idx_local] = min(dists);
nearest_idx = search_range(idx_local);

% 路径点平均间距
if size(path,1) > 1
    avg_spacing = mean(sqrt(sum(diff(path(:,1:2)).^2, 2)));
else
    avg_spacing = 0.5;
end

% 从最近点向前看 Delta
lookahead_steps = max(1, round(Delta / avg_spacing));
target_idx = min(nearest_idx + lookahead_steps, size(path,1));
target_idx = max(target_idx, 2);

Pk = path(target_idx, :)';

% 期望航向
dx = Pk(1) - x;
dy = Pk(2) - y;
psi_raw = atan2(dy, dx);

% 侧滑角补偿
if abs(u) > 0.1
    beta = atan2(v, u);
else
    beta = 0;
end
psi_raw = psi_raw - beta;

% 航向连续化
dpsi = psi_raw - psi_continuous;
while dpsi >  pi, dpsi = dpsi - 2*pi; end
while dpsi < -pi, dpsi = dpsi + 2*pi; end
psi_continuous = psi_continuous + dpsi;

psi_d = mod(psi_continuous + pi, 2*pi) - pi;
end