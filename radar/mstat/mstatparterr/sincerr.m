clc; clear; close all;

%% 1. 参数设置
M = 64; 
N = 32;
h_amp = 1.66135;       % 幅值
tau_true = 0.76330;    % 真实时延
tau_est  = 0.8;    % 估计时延 (仅差 0.00001)
nu_idx = 0;            % 假设多普勒为0 (只观察时延切片)

% 设置观察窗口 (以目标为中心，前后几个点)
plot_center = 1;       % 目标大约在 index 1 (0.76)
span = 5;              % 左右显示的范围
x_range = (plot_center - span) : 0.1 : (plot_center + span);

%% 2. 生成离散点 (Simulation Points)
% 我们构造一个 dd 为单位冲激，这样 generate_echo 输出的就是信道响应本身
dd = zeros(M, N);
dd(1,1) = 1; 

% 调用你的函数生成 M x N 矩阵
Y_true_mat = generate_echo_physical_truth(dd, h_amp, tau_true, nu_idx, M, N);
Y_est_mat  = generate_echo_physical_truth(dd, h_amp, tau_est,  nu_idx, M, N);

% 提取第 1 列 (对应 Doppler = 0)
y_true_discrete = Y_true_mat(:, 1);
y_est_discrete  = Y_est_mat(:, 1);
y_resid_discrete = y_true_discrete - y_est_discrete;

%% 3. 生成平滑曲线 (Analytic Curves) 用于绘图美观
% 直接利用函数内的 Dirichlet Kernel 公式计算高分辨率曲线
L_fine = plot_center - span : 0.02 : plot_center + span; % 细网格

% 定义计算函数 (从你的代码中提取的核心数学逻辑)
calc_response = @(l_axis, tau) arrayfun(@(l) compute_point_val(l, tau, M, h_amp), l_axis);

curve_true = calc_response(L_fine, tau_true);
curve_est  = calc_response(L_fine, tau_est);
curve_resid = curve_true - curve_est;

%% 4. 绘图 (仿照参考图风格)
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; grid on; box on;

% 绘制平滑曲线 (实部)
p1 = plot(L_fine, real(curve_true), 'b--', 'LineWidth', 2, 'DisplayName', '真实目标 (Exact)');
p2 = plot(L_fine, real(curve_est),  'g-.', 'LineWidth', 2, 'DisplayName', '估计目标 (Inexact)');
p3 = plot(L_fine, real(curve_resid),'r-',  'LineWidth', 2, 'DisplayName', '相减残差 (Residual)');

% 绘制离散采样点 (你的系统实际看到的点)
% 注意：MATLAB 索引是 1~M，而物理公式通常 0~M-1。
% 你的函数 generate_echo 内部 ndgrid(0:M-1)，所以对应 MATLAB 的 1:M
indices = 1:M; 
% 只画视野范围内的点
mask = indices >= min(x_range) & indices <= max(x_range);

scatter(indices(mask)-1, real(y_true_discrete(mask)), 50, 'b', 'filled', 'HandleVisibility', 'off');
scatter(indices(mask)-1, real(y_est_discrete(mask)),  50, 'g', 'filled', 'HandleVisibility', 'off');
scatter(indices(mask)-1, real(y_resid_discrete(mask)),50, 'r', 'filled', 'HandleVisibility', 'off');

% 标注与美化
xlabel('Delay Index (0-based)');
ylabel('Amplitude (Real Part)');
title(sprintf('Sinc Subtraction Analysis\nTrue: %.5f | Est: %.5f | Diff: %.1e,绘制的是实部', ...
    tau_true, tau_est, abs(tau_true - tau_est)));
legend('Location', 'best');

% 调整坐标轴
xlim([min(L_fine), max(L_fine)]);
ylim([-0.5, h_amp * 1.2]); 

% 添加一条零线
yline(0, 'k-', 'HandleVisibility', 'off', 'Alpha', 0.2);

hold off;

%% --- 辅助函数：计算单点物理值 (核心数学逻辑) ---
function val = compute_point_val(l_grid_val, tau, M, h)
    % 模拟你的代码中的 l_val 处理
    l_val = l_grid_val - tau;
    
    % 【关键】周期性 Wrap 处理 (你的代码中的 round 逻辑)
    l_val = l_val - M * round(l_val / M);
    
    % Dirichlet Kernel 计算
    if abs(sin(pi * l_val / M)) < 1e-9
        F_term = 1;
    else
        % 你的公式: exp(...) * (sin / (M*sin))
        F_term = exp(-1i*pi*(M-1)*l_val/M) .* (sin(pi*l_val)./(M*sin(pi*l_val/M)));
    end
    
    val = h * F_term;
end

%% --- 你的原始函数 (作为本地函数被调用) ---
function Y_out = generate_echo_physical_truth(dd, h, tau_idx, nu_idx, M, N)
    % Dirichlet Kernel 物理真值生成
    [L_grid, K_grid] = ndgrid(0:M-1, 0:N-1);
    
    % Doppler Kernel (G_term)
    k_val = K_grid - nu_idx;
    k_val = k_val - N * round(k_val / N);
    mask_peak = abs(sin(pi * k_val / N)) < 1e-9;
    G_term = zeros(size(k_val));
    G_term(mask_peak) = 1;
    mn = ~mask_peak;
    if any(mn(:))
        kv = k_val(mn);
        G_term(mn) = exp(-1i*pi*(N-1)*kv/N) .* (sin(pi*kv)./(N*sin(pi*kv/N)));
    end
    
    % Delay Kernel (F_term)
    l_val = L_grid - tau_idx;
    l_val = l_val - M * round(l_val / M);
    mask_peak = abs(sin(pi * l_val / M)) < 1e-9;
    F_term = zeros(size(l_val));
    F_term(mask_peak) = 1;
    mn = ~mask_peak;
    if any(mn(:))
        lv = l_val(mn);
        F_term(mn) = exp(-1i*pi*(M-1)*lv/M) .* (sin(pi*lv)./(M*sin(pi*lv/M)));
    end
    
    theta_term = exp(-1j * 2 * pi * tau_idx * nu_idx / (M*N));
    hw = h * G_term .* F_term * theta_term;
    
    % 2D 循环卷积
    % 如果 dd 是单位冲激，这里输出的其实就是 hw
    Y_out = ifft2(fft2(hw) .* fft2(dd));
end