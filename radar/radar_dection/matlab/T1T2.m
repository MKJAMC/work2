% 生成：最大/最小增益比 < 5，且 T1 和 T2 时延多普勒间隔均 < 1
clc; clear; close all;
%% ================= 参数设置 =================
% 基础参数 (保持 resolution.m 与仿真设置一致)
M = 64;
N = 32;
fc = 64e9;
delta_f = 120e3;
c = 3e8;
delay_res = 1 / (M * delta_f);
doppler_res = 1 / (N * (1/delta_f));
distance_res = c * delay_res / 2;
v_res = doppler_res * c / 2 / fc;

% OTFS 调制参数
l_max = 20; k_max = 2;
cen_modu = 2; guard_modu = 4;
bits_per_qam_symbol = log2(cen_modu);
bits_per_guard_symbol = log2(guard_modu);
g = 4;

% 通信距离通常是雷达系统最大检测距离的两倍，通信系统只需要单程传播
R_max = (l_max/2) * distance_res;
V_max = N * v_res / 2;

% 仿真配置
num_targets = 3;
SNR_vec = 10:2:18;  % SNR 范围 (dB)
%% SNR
num_monte_carlo = 200; % 蒙特卡洛次数 (为了画出平滑曲线设为 200)

% 预计算索引参数
guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
data_rows = (l_max + 2):(M - l_max);
guard_rows = [1:(l_max + 1), (M - l_max + 1):M];

% 结果存储
MSE_Range = zeros(length(SNR_vec), 1);
MSE_Velocity = zeros(length(SNR_vec), 1);
Detection_Prob = zeros(length(SNR_vec), 1); % 新增检测概率统计

fprintf('======================================================\n');
fprintf('开始 Monte Carlo 仿真 (T1/T2 不可分辨版 + 增益比<5)\n');
fprintf('SNR 范围: [%s] dB\n', num2str(SNR_vec));
fprintf('======================================================\n');

%% ================= 主循环 (SNR) =================
total_timer = tic;

for s_idx = 1:length(SNR_vec)
    current_SNR = SNR_vec(s_idx);
    fprintf('\n正在仿真 SNR = %d dB ...\n', current_SNR);

    % --- 功率计算逻辑 ---
    power_persym = 10.^(current_SNR/10);
    ans2 = (M - 2*l_max - 1)*N + (N - 4*k_max - 1)*(2*l_max + 1);
    power = power_persym * ans2;
    ans3 = 430;
    factor0 = sqrt(power ./ ans3);

    % 累加误差变量
    total_valid_targets_SNR = 0; % 关键：每个 SNR 点开始前，有效目标计数清零
    sum_sq_err_range = 0;
    sum_sq_err_vel = 0;
    dist = zeros(1, num_monte_carlo);

    %% 蒙特卡洛循环
    for mc = 1:num_monte_carlo
        % mc=26;
        disp(mc);
        
        rng(mc); % 确保复现

        % ================= 【核心修改区：对齐 nonscatter_test.m】 =================
        target_amp_ratio_min = 1.2;  % 改为 1.2
        target_amp_ratio_max = 5.0;  

        % === 定义间隔约束 ===
        % 1. T1 和 T2 的距离、速度间隔均 < 1 倍分辨率 (取 0.3~0.6 之间)
        gap_R12_min = 0.3 * distance_res;
        gap_R12_max = 0.6 * distance_res;
        gap_V12_min = 0.3 * v_res;
        gap_V12_max = 0.6 * v_res;

        % 2. T3 与前两个目标的最小间隔 > 1.2 倍分辨率 (保证 T3 是独立分散的)
        safe_gap_R3  = 1.2 * distance_res;
        safe_gap_V3  = 1.2 * v_res;

        is_valid = false;
        loop_safety_count = 0;

        while ~is_valid
            loop_safety_count = loop_safety_count + 1;

            % --- A. 生成距离 (Range) ---
            R1_min = 2.5 * distance_res;
            R1_max = R_max * 0.35;
            R1 = R1_min + rand * (R1_max - R1_min);

            % 生成 R2 (与 R1 间隔 < 1)
            actual_gap_R12 = gap_R12_min + rand * (gap_R12_max - gap_R12_min);
            R2 = R1 + actual_gap_R12;

            % --- B. 【定向生成 R3】 ---
            desired_ratio = target_amp_ratio_min + rand * (target_amp_ratio_max - target_amp_ratio_min);
            R3 = R1 * sqrt(desired_ratio);

            % 验证 R3 是否满足安全间隔和最大检测距离约束
            if (R3 > R2 + safe_gap_R3) && (R3 < R_max * 0.95)
                T_Range = [R1, R2, R3];

                % --- C. 生成速度 (Velocity) ---
                V1 = (rand * 2 - 1) * V_max * 0.7;

                % 生成 V2
                actual_gap_V12 = gap_V12_min + rand * (gap_V12_max - gap_V12_min);
                if rand > 0.5
                    V2 = V1 + actual_gap_V12;
                else
                    V2 = V1 - actual_gap_V12;
                end

                valid_vel3 = false;
                vel_retry = 0;
                while ~valid_vel3 && vel_retry < 50
                    V3 = (rand * 2 - 1) * V_max * 0.9;
                    if abs(V3 - V1) > safe_gap_V3 && abs(V3 - V2) > safe_gap_V3
                        valid_vel3 = true;
                        T_Velocity = [V1, V2, V3];
                    end
                    vel_retry = vel_retry + 1;
                end

                if valid_vel3
                    % --- D. 计算幅值并最后确认 ---
                    raw_trend = 1 ./ (T_Range.^2);
                    current_power = mean(raw_trend.^2);
                    true_h_abs = raw_trend / sqrt(current_power);

                    current_dyn_range = max(true_h_abs) / min(true_h_abs);

                    % 【最关键的修改：与 nonscatter_test.m 的逻辑完全一致】
                    % 去掉了 >= target_amp_ratio_min 的判断，只保留 <= target_amp_ratio_max
                    if current_dyn_range <= target_amp_ratio_max
                        is_valid = true;
                        true_h = true_h_abs .* exp(1j * 2 * pi * rand(1, num_targets));
                    end
                end
            end

            if loop_safety_count > 5000
                R1_min = R1_min * 1.05; % 死循环时放宽下限
                loop_safety_count = 0;
            end
        end
        % ===================================================

        true_li = T_Range ./ distance_res;
        true_ki = T_Velocity ./ v_res;

        % --- 2. 生成发射信号 dd ---
        dd = zeros(M, N);
        dd(1,1) = 1000; % 导频

        % 数据区
        num_groups_data = N / g;
        for row = data_rows
            for j = 1:num_groups_data
                im_bits = randi([0 1], 1, log2(nchoosek(g, 1)));
                qam_bits = randi([0 1], 1, bits_per_qam_symbol);
                active_idx = bi2de(im_bits, 'left-msb') + 1;
                sym_val = qammod(bi2de(qam_bits, 'left-msb'), cen_modu, 'UnitAveragePower', true);
                dd(row, (j-1)*g + active_idx) = sym_val * factor0;
            end
        end
        % 保护区
        num_guard_cols = length(guard_cols_range);
        for row = guard_rows
            i_col = 1;
            while i_col <= num_guard_cols
                grp_sz = min(g, num_guard_cols - i_col + 1);
                if grp_sz > 1
                    if grp_sz == g
                        act_idx = bi2de(randi([0 1], 1, floor(log2(nchoosek(g, 1)))), 'left-msb') + 1;
                    else
                        act_idx = randi(grp_sz);
                    end
                    sym_val = qammod(bi2de(randi([0 1], 1, bits_per_guard_symbol), 'left-msb'), guard_modu, 'UnitAveragePower', true);
                    dd(row, guard_cols_range(i_col + act_idx - 1)) = sym_val * factor0;
                end
                i_col = i_col + grp_sz;
            end
        end

        % --- 3. 生成接收信号 (Ground Truth) ---
        Y_clean = zeros(M, N);
        for p = 1:num_targets
            Y_clean = Y_clean + generate_echo_physical_truth(dd, true_h(p), true_li(p), true_ki(p), M, N);
        end
        noise = sqrt(1/2) * (randn(M, N) + 1j * randn(M, N));
        Y = Y_clean + noise;

        % --- 4. MTASA 检测算法 ---
        [est_li, est_ki] = run_mtasa_detection(Y, dd, M, N, num_targets, delta_f);

        % --- 5. 误差计算 ---
        est_dist = est_li * distance_res;
        est_vel  = est_ki * v_res;

        % 1. 构建误差矩阵
        dist_err_mat = zeros(num_targets, num_targets);
        vel_err_mat  = zeros(num_targets, num_targets);

        for t = 1:num_targets
            for e = 1:num_targets
                % 距离误差 (考虑 M 周期性)
                diff_l = est_li(e) - true_li(t);
                diff_l = diff_l - M * round(diff_l / M);
                dist_err_mat(t, e) = abs(diff_l * distance_res);

                % 速度误差 (考虑 N 周期性)
                diff_k = est_ki(e) - true_ki(t);
                diff_k = diff_k - N * round(diff_k / N);
                vel_err_mat(t, e) = abs(diff_k * v_res);
            end
        end

        % 2. 贪心匹配
        cost_mat = dist_err_mat;
        matched_pairs = [];

        for pair_idx = 1:num_targets % 强制匹配所有目标
            [min_val, min_linear_idx] = min(cost_mat(:));
            [row, col] = ind2sub(size(cost_mat), min_linear_idx);
            matched_pairs = [matched_pairs; row, col];
            cost_mat(row, :) = inf;
            cost_mat(:, col) = inf;
        end

        % 3. 统计误差与检测成功数
        mc_range_sq_err = 0;
        mc_vel_sq_err = 0;
   

        % 定义检测成功的判定门限（仅用于计算检测概率，不影响 MSE）
        detect_thresh_dist = 1.0 * distance_res;
        detect_thresh_vel  = 1.0 * v_res;
        is_sample_correct = true; % 默认当前样本完全正确
        
        for i = 1:num_targets
            t_idx = matched_pairs(i, 1);
            e_idx = matched_pairs(i, 2);
            r_err = dist_err_mat(t_idx, e_idx);
            v_err = vel_err_mat(t_idx, e_idx);
            
            % 累加误差 (供MSE计算使用)
            mc_range_sq_err = mc_range_sq_err + r_err^2;
            mc_vel_sq_err   = mc_vel_sq_err + v_err^2;
            
            % 【修改核心】只要有1个目标的距离或速度误差超出门限，整个样本判定为错
            if r_err >= detect_thresh_dist || v_err >= detect_thresh_vel
                is_sample_correct = false;
            end
        end
        
        total_valid_targets_SNR = total_valid_targets_SNR + num_targets;
        sum_sq_err_range = sum_sq_err_range + mc_range_sq_err;
        sum_sq_err_vel   = sum_sq_err_vel + mc_vel_sq_err;
        
        if ~exist('total_success_count_SNR', 'var') || mc == 1
            total_success_count_SNR = 0;
        end
        
        % 如果该样本3个目标全都找对了，正确样本数才 + 1
        if is_sample_correct
            total_success_count_SNR = total_success_count_SNR + 1;
        end


        dist(mc) = sqrt(mc_range_sq_err);



        %% 打印某个蒙特卡罗下的具体数据
        % num_targets_count = length(true_li);
        % % 1. 打印距离索引信息
        % fprintf('\n--- 距离维度 (Range Index) ---\n');
        % % 适当放宽幅值列的宽度以容纳复数格式
        % fprintf('%-10s | %-12s | %-12s | %-10s | %-20s\n', '目标', '真实值', '估计值', '偏差', '幅值 (复数)');
        % for i = 1:num_targets_count
        %     % 匹配对应的估计索引
        %     fprintf('Path [%d]   | %12.5f | %12.5f | %+10.5f | %8.5f%+8.5fi\n', ...
        %         i, true_li(i), est_li(i), est_li(i) - true_li(i), real(true_h(i)), imag(true_h(i)));
        % end

        % % 2. 打印速度维度信息
        % fprintf('\n--- 速度维度 (Velocity Index) ---\n');
        % fprintf('%-10s | %-12s | %-12s | %-10s | %-20s\n', '目标', '真实值', '估计值', '偏差', '幅值 (复数)');
        % for i = 1:num_targets_count
        %     fprintf('Path [%d]   | %12.5f | %12.5f | %+10.5f | %8.5f%+8.5fi\n', ...
        %         i, true_ki(i), est_ki(i), est_ki(i) - true_ki(i), real(true_h(i)), imag(true_h(i)));
        % end  

        % figure('Name', 'DD Domain Range-Doppler Heatmap', 'Color', 'w');

        % % 1. 绘制接收信号幅度热力图
        % % X轴: Doppler (1~N), Y轴: Delay (1~M)
        % imagesc(1:N, 1:M, abs(Y));
        % colormap('jet'); % 使用高对比度颜色
        % colorbar;
        % xlabel('Doppler Index (Column Index)');
        % ylabel('Delay Index (Row I ndex)');
        % title(['Received Signal Y (Magnitude) & Ground Truths (Red X) | SNR=' num2str(current_SNR) 'dB']);
        % hold on;

        % % 2. 循环绘制真实目标位置
        % for t = 1:num_targets
        %     % --- 坐标转换逻辑 ---

        %     % A. 距离坐标 (Delay)
        %     % 真实值是 0 ~ M-1，MATLAB索引是 1 ~ M
        %     y_pos = true_li(t) + 1;

        %     % B. 速度坐标 (Doppler) - 关键！
        %     % 物理真值可能是负数 (例如 -3.5)
        %     % 在矩阵存储中，负频率折叠到了后半部分 (N/2 ~ N)
        %     k_val = true_ki(t);

        %     if k_val < 0
        %         % 如果是负速度，映射到 N + k_val + 1
        %         % 例如 N=32, k=-2 => index = 32 - 2 + 1 = 31
        %         x_pos = N + k_val + 1;
        %     else
        %         % 如果是正速度，直接 + 1
        %         % 例如 k=2 => index = 2 + 1 = 3
        %         x_pos = k_val + 1;
        %     end

        %     % 3. 绘制红色的 'x'
        %     plot(x_pos, y_pos, 'rx', 'MarkerSize', 12, 'LineWidth', 2.5);

        %     % (可选) 添加目标编号标签
        %     text(x_pos, y_pos - 1, sprintf(' T%d', t), ...
        %         'Color', 'white', 'FontWeight', 'bold', 'FontSize', 10);
        % end
        % hold off;
        % % 强制让Y轴方向符合矩阵直观 (行号从上到下增加)
        % set(gca, 'YDir', 'reverse');
        % % % disp(sqrt(dist(mc)));


    end

    if total_valid_targets_SNR > 0
        MSE_Range(s_idx) = sqrt(sum_sq_err_range / total_valid_targets_SNR);
        MSE_Velocity(s_idx) = sqrt(sum_sq_err_vel / total_valid_targets_SNR);
    else
        MSE_Range(s_idx) = NaN;
        MSE_Velocity(s_idx) = NaN;
    end

    Detection_Prob(s_idx) = total_success_count_SNR / num_monte_carlo;

    fprintf(' 完成! 有效样本数: %d, MSE_Range: %.4e, 检测概率: %.2f%%\n', ...
        total_valid_targets_SNR, MSE_Range(s_idx), Detection_Prob(s_idx)*100);

    figure; %检查在某个snr下所有蒙特卡洛仿真出的距离误差
    plot(1:num_monte_carlo, dist, '-o', 'LineWidth', 1.5);
    grid on;
    xlabel('Monte Carlo Iteration Index');
    ylabel('Dist Range (m^2)');
    title(['Per-Iterat  ion Range Error at SNR = ' num2str(current_SNR) ' dB']);
end
toc(total_timer);

%% ================= 绘图 =================
figure('Color', 'w', 'Position', [100, 100, 1200, 400]);

subplot(1,3,1);
semilogy(SNR_vec, MSE_Range, '-ro', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)'); ylabel('RMSE of Range (m)');
title('Range Estimation Accuracy (RMSE)');

subplot(1,3,2);
semilogy(SNR_vec, MSE_Velocity, '-bs', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)'); ylabel('RMSE of Velocity (m/s)');
title('Velocity Estimation Accuracy (RMSE)');

subplot(1,3,3);
plot(SNR_vec, Detection_Prob * 100, '-kP', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)'); ylabel('Detection Probability (%)');
title('Target Detection Performance');
ylim([0 105]);

%% ================= 核心算法函数 =================
function [est_li, est_ki] = run_mtasa_detection(Y, dd, M, N, num_targets, delta_f)
est_li = zeros(1, num_targets);
est_ki = zeros(1, num_targets);
est_h  = zeros(1, num_targets);

Residual = Y;
for k = 1:num_targets
    [h_hat, l_hat, k_hat] = two_stage_search_matched_golden(Residual, dd, M, N);
    est_li(k) = l_hat;
    est_ki(k) = k_hat;
    est_h(k)  = h_hat;
    Signal_Est = generate_echo_physical_truth(dd, h_hat, l_hat, k_hat, M, N);
    Residual = Residual - Signal_Est;
end

MaxIter = 10;
for iter = 1:MaxIter
    prev_li = est_li;
    prev_ki = est_ki;

    for k = 1:num_targets
        Interference = zeros(M, N);
        for j = 1:num_targets
            if j ~= k
                Interference = Interference + generate_echo_physical_truth(dd, est_h(j), est_li(j), est_ki(j), M, N);
            end
        end

        Y_clean_k = Y - Interference;
        [h_new, l_new, k_new] = two_stage_search_matched_golden(Y_clean_k, dd, M, N);

        est_li(k) = l_new;
        est_ki(k) = k_new;
        est_h(k)  = h_new;
    end

    if norm([est_li - prev_li, est_ki - prev_ki]) < 1e-3
        break;
    end
end

mask_neg = est_ki > (N/2);
est_ki(mask_neg) = est_ki(mask_neg) - N;
est_li = mod(est_li, M);
end

function [h_hat, l_hat, k_hat] = two_stage_search_matched_golden(Y_received, dd, M, N)
max_corr = -1; l_coarse_idx = 0; k_coarse_idx = 0;

for l = 0:M-1
    for k = 0:N-1
        dd_shifted = circshift(dd, [l, k]);
        corr_val = abs(sum(sum(Y_received .* conj(dd_shifted))))^2;
        if corr_val > max_corr
            max_corr = corr_val;
            l_coarse_idx = l; k_coarse_idx = k;
        end
    end
end

if k_coarse_idx > N/2
    k_coarse_idx = k_coarse_idx - N;
end

delta_f = 120e3;
T = 1 / delta_f;

tau_c = l_coarse_idx / (M * delta_f);
nu_c  = k_coarse_idx / (N * T);

res_tau = 1 / (M * delta_f);
res_nu  = 1 / (N * T);

a_l = tau_c - res_tau; a_u = tau_c + res_tau;
b_l = nu_c - res_nu;   b_u = nu_c + res_nu;

mu = (sqrt(5) - 1) / 2;
Iter = 10;

for i = 1:Iter
    I_a = a_u - a_l;
    I_b = b_u - b_l;

    a1 = a_l + (1 - mu) * I_a; a2 = a_l + mu * I_a;
    b1 = b_l + (1 - mu) * I_b; b2 = b_l + mu * I_b;

    vals = zeros(2, 2);
    E11 = generate_echo_physical_truth_continuous(dd, 1, a1, b1, M, N);
    vals(1,1) = abs(sum(sum(Y_received .* conj(E11))))^2;
    E12 = generate_echo_physical_truth_continuous(dd, 1, a1, b2, M, N);
    vals(1,2) = abs(sum(sum(Y_received .* conj(E12))))^2;
    E21 = generate_echo_physical_truth_continuous(dd, 1, a2, b1, M, N);
    vals(2,1) = abs(sum(sum(Y_received .* conj(E21))))^2;
    E22 = generate_echo_physical_truth_continuous(dd, 1, a2, b2, M, N);
    vals(2,2) = abs(sum(sum(Y_received .* conj(E22))))^2;

    [~, max_idx] = max(vals(:));
    [r_idx, c_idx] = ind2sub([2, 2], max_idx);

    if r_idx == 1 && c_idx == 1
        a_u = a2; b_u = b2;
    elseif r_idx == 1 && c_idx == 2
        a_u = a2; b_l = b1;
    elseif r_idx == 2 && c_idx == 1
        a_l = a1; b_u = b2;
    elseif r_idx == 2 && c_idx == 2
        a_l = a1; b_l = b1;
    end
end

tau_hat_final = (a_l + a_u) / 2;
nu_hat_final  = (b_l + b_u) / 2;

l_hat = tau_hat_final * (M * delta_f);
k_hat = nu_hat_final * (N * T);

E_final = generate_echo_physical_truth_continuous(dd, 1, tau_hat_final, nu_hat_final, M, N);
h_hat = sum(sum(Y_received .* conj(E_final))) / sum(sum(abs(E_final).^2));
end

function Y_out = generate_echo_physical_truth_continuous(dd, h, tau_phy, nu_phy, M, N)
delta_f = 120e3;
T = 1/delta_f;
tau_idx = tau_phy * (M * delta_f);
nu_idx  = nu_phy * (N * T);
Y_out = generate_echo_physical_truth(dd, h, tau_idx, nu_idx, M, N);
end

function Y_out = generate_echo_physical_truth(dd, h, tau_idx, nu_idx, M, N)
[L_grid, K_grid] = ndgrid(0:M-1, 0:N-1);
k_val = K_grid - nu_idx;
mask_peak = abs(sin(pi * k_val / N)) < 1e-9;
G_term = zeros(size(k_val));
G_term(mask_peak) = 1;
mn = ~mask_peak;
if any(mn(:))
    kv = k_val(mn);
    G_term(mn) = exp(-1i*pi*(N-1)*kv/N) .* (sin(pi*kv)./(N*sin(pi*kv/N)));
end
l_val = L_grid - tau_idx;
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
Y_out = ifft2(fft2(hw) .* fft2(dd));
end