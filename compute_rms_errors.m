function [rms_x, rms_y, rms_3d, max_err] = compute_rms_errors(log_eta, path)
% 计算三方向 RMS 误差与三维合成 RMS
N = size(log_eta, 2);
errors_x = zeros(1, N);
errors_y = zeros(1, N);
errors_3d = zeros(1, N);
for k = 1:N
    [~, idx] = min(abs(path(:,1) - log_eta(1,k)));
    errors_x(k) = abs(path(idx,1) - log_eta(1,k));
    errors_y(k) = abs(path(idx,2) - log_eta(2,k));
    errors_3d(k) = sqrt(errors_x(k)^2 + errors_y(k)^2);
end
rms_x = sqrt(mean(errors_x.^2));
rms_y = sqrt(mean(errors_y.^2));
rms_3d = sqrt(mean(errors_3d.^2));
max_err = max(errors_3d);
fprintf('==========================================\n');
fprintf('X 方向 RMS: %.4f m\n', rms_x);
fprintf('Y 方向 RMS: %.4f m\n', rms_y);
fprintf('二维合成 RMS: %.4f m\n', rms_3d);
fprintf('最大误差: %.4f m\n', max_err);
fprintf('==========================================\n');
end