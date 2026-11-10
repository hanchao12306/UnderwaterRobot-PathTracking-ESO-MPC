function features = analyze_path(path)
% ============================================================
% 稳健版路径特征分析
% - 先去除重复点
% - 用多尺度平滑计算曲率
% ============================================================

% ---------- 步骤1：去除重复点 ----------
min_dist = 0.05;   % 最小点间距 5cm
keep_idx = [1];
for i = 2:size(path,1)
    if norm(path(i,:) - path(keep_idx(end),:)) > min_dist
        keep_idx = [keep_idx, i];
    end
end
path_clean = path(keep_idx, :);
fprintf('原始路径点数: %d, 去重后: %d\n', size(path,1), size(path_clean,1));

% ---------- 步骤2：均值滤波（进一步平滑）----------
window = 5;
path_smooth = path_clean;
for j = 1:3
    path_smooth(:,j) = smoothdata(path_clean(:,j), 'gaussian', window);
end

% ---------- 步骤3：计算切向量 ----------
N = size(path_smooth, 1);
tangent = zeros(N-1, 3);
seg_len = zeros(N-1, 1);
for i = 1:N-1
    d = path_smooth(i+1,:) - path_smooth(i,:);
    seg_len(i) = norm(d);
    if seg_len(i) > 1e-6
        tangent(i,:) = d / seg_len(i);
    else
        tangent(i,:) = [1, 0, 0];   % 防零向量
    end
end

% ---------- 步骤4：多尺度曲率计算 ----------
% 用间隔 5 个点的差分，避免局部噪声
stride = 5;
curvature = zeros(N-2*stride, 1);
for i = 1:N-2*stride
    T1 = path_smooth(i+stride,:) - path_smooth(i,:);
    T2 = path_smooth(i+2*stride,:) - path_smooth(i+stride,:);
    if norm(T1) > 1e-6 && norm(T2) > 1e-6
        T1 = T1 / norm(T1);
        T2 = T2 / norm(T2);
        dtheta = acos(max(min(dot(T1,T2), 1), -1));
        ds = (norm(path_smooth(i+stride,:) - path_smooth(i,:)) + ...
              norm(path_smooth(i+2*stride,:) - path_smooth(i+stride,:))) / 2;
        curvature(i) = dtheta / ds;
    end
end

% ---------- 步骤5：滤除异常曲率 ----------
kappa_valid = curvature(curvature > 1e-6 & curvature < 10);   % 合理范围
if isempty(kappa_valid)
    kappa_valid = 0.05;   % 默认
end

% ---------- 步骤6：路径特征统计 ----------
L_total = sum(seg_len);
kappa_max = max(kappa_valid);
kappa_mean = mean(kappa_valid);
kappa_median = median(kappa_valid);
R_min = 1 / (kappa_max + 1e-6);       % 最小曲率半径
R_median = 1 / (kappa_median + 1e-6);  % 中位曲率半径
R_mean = 1 / (kappa_mean + 1e-6);

% ---------- 步骤7：路径分段分析（每 5 米一段）----------
num_seg = ceil(L_total / 5);
seg_R_min = zeros(num_seg, 1);
seg_start = 1;
cum_len = 0;
for s = 1:num_seg
    cum_len = cum_len + 5;
    idx_end = find(cumsum(seg_len) > cum_len, 1);
    if isempty(idx_end), idx_end = N-1; end
    if idx_end > seg_start + 2
        seg_kappa = curvature(seg_start:min(idx_end, length(curvature)));
        seg_kappa = seg_kappa(seg_kappa > 1e-6 & seg_kappa < 10);
        if ~isempty(seg_kappa)
            seg_R_min(s) = 1 / max(seg_kappa);
        else
            seg_R_min(s) = 10;   % 默认直线
        end
    end
    seg_start = idx_end;
end
seg_R_min = seg_R_min(seg_R_min > 0);

% ---------- 步骤8：Z方向斜率 ----------
slope_z = zeros(N-1, 1);
for i = 1:N-1
    if seg_len(i) > 1e-6
        slope_z(i) = abs(path_smooth(i+1,3) - path_smooth(i,3)) / seg_len(i);
    end
end
max_slope_z = max(slope_z);
mean_slope_z = mean(slope_z);

% ---------- 步骤9：输出 ----------
features.L_total = L_total;
features.R_min = R_min;
features.R_median = R_median;
features.R_mean = R_mean;
features.kappa_max = kappa_max;
features.kappa_mean = kappa_mean;
features.max_slope_z = max_slope_z;
features.mean_slope_z = mean_slope_z;
features.seg_R_min = seg_R_min;
features.path_clean = path_smooth;

% ---------- 步骤10：打印 ----------
fprintf('========= 路径特征分析（稳健版）=========\n');
fprintf('总弧长 L          : %.2f m\n', L_total);
fprintf('最小曲率半径 R_min: %.3f m\n', R_min);
fprintf('中位曲率半径       : %.3f m\n', R_median);
fprintf('平均曲率半径       : %.3f m\n', R_mean);
fprintf('最大曲率 kappa_max: %.4f 1/m\n', kappa_max);
fprintf('平均曲率 kappa_mean: %.4f 1/m\n', kappa_mean);
fprintf('Z方向最大斜率      : %.3f\n', max_slope_z);
fprintf('Z方向平均斜率      : %.3f\n', mean_slope_z);
fprintf('各 5 米段 R_min:\n');
for s = 1:length(seg_R_min)
    fprintf('  段%d: R_min = %.2f m\n', s, seg_R_min(s));
end
fprintf('=========================================\n');
end