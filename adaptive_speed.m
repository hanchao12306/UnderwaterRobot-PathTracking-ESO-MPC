function [ud, wd] = adaptive_speed(eta, path, target_idx, params)
% ============================================================
% 自适应纵荡速度 + 垂向制导（论文4.2.3、4.2.4节）
% 【关键修复】用二维最近点（x,y）代替三维最近点
% 避免 z 方向的偏离干扰最近点的选择
% ============================================================

x = eta(1); y = eta(2); z = eta(3);

% 【修复】只用 x,y 计算最近点，忽略 z
dists_xy = sqrt((path(:,1) - x).^2 + (path(:,2) - y).^2);
[~, nearest_idx] = min(dists_xy);
Pk = path(nearest_idx, :)';

% 横向偏差（式4-12）
e_cross = sqrt((x - Pk(1))^2 + (y - Pk(2))^2);
% 垂向偏差（式4-13）
e_z = abs(z - Pk(3));

% 自适应纵荡速度（式4-15）
ud = params.umax - (params.umax - params.umin) * ...
     tanh(params.ke * sqrt(e_cross^2 + e_z^2));

% 垂向制导（式4-10），用二维最近点的 z
wd = params.kz * (Pk(3) - z);
wd = max(min(wd, params.wmax), -params.wmax);
end