% 对第三个目标随机生成，保证幅值比>50即可。
%% 644s
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
SNR_vec = 10:2:18 ;  % SNR 范围 (dB)
%% SNR
num_monte_carlo = 200; % 蒙特卡洛次数


% 预计算索引参数
guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
data_rows = (l_max + 2):(M - l_max);
guard_rows = [1:(l_max + 1), (M - l_max + 1):M];

% 结果存储
MSE_Range = zeros(length(SNR_vec), 1);
MSE_Velocity = zeros(length(SNR_vec), 1);

fprintf('======================================================\n');
fprintf('开始 Monte Carlo 仿真 (终极修正版：黄金分割法 + 蒙特卡洛修复)\n');
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
    mc_sq_err_history = zeros(1, num_monte_carlo);
    dist=[];
    %% 蒙特卡洛循环
    for mc = 1:num_monte_carlo
        disp(mc);
        % mc=443;
        % --- 1. 生成随机目标 (与 scatter_test.m 保持绝对一致) ---
        rng(mc);
        max_amp_ratio = 5;       % 最大/最小回波增益比 < 5
        
        % === 定义间隔约束 (针对所有目标) ===
        gap_R_min = 1.2 * distance_res; % 距离间隔 > 1 倍分辨率 (设1.2留容差)
        gap_V_min = 1.2 * v_res;        % 速度间隔 > 1 倍分辨率
        
        is_valid = false;
        loop_safety_count = 0;
        
        while ~is_valid
            loop_safety_count = loop_safety_count + 1;
            
            % === A. 生成距离 (Range) ===
            % R1 必须够远，否则塞不下 R2 和 R3，同时满足 R3/R1 < sqrt(5)
            R1_min = 2.5 * distance_res; 
            R1_max = R_max * 0.4; 
            R1 = R1_min + rand * (R1_max - R1_min);
            
            % R2 与 R1 保持间隔
            R2_min = R1 + gap_R_min;
            R2_max = R2_min + 1.0 * distance_res;
            R2 = R2_min + rand * (R2_max - R2_min);
            
            % R3 与 R2 保持间隔，且满足增益上限约束
            R3_min = R2 + gap_R_min;
            R3_max_dyn = R1 * sqrt(max_amp_ratio) * 0.98; 
            
            if R3_max_dyn > R3_min
                R3_max = min(R_max * 0.95, R3_max_dyn);
                R3 = R3_min + rand * (R3_max - R3_min);
                T_Range = [R1, R2, R3];
                
                % === B. 生成速度 (Velocity) ===
                valid_vel = false;
                vel_retry = 0;
                while ~valid_vel && vel_retry < 50
                    % 三个目标速度全随机
                    V1 = (rand * 2 - 1) * V_max * 0.8;
                    V2 = (rand * 2 - 1) * V_max * 0.8;
                    V3 = (rand * 2 - 1) * V_max * 0.8;
                    
                    % 必须满足任意两两之间的速度间隔 > 1倍多普勒分辨率
                    if abs(V1 - V2) > gap_V_min && abs(V2 - V3) > gap_V_min && abs(V1 - V3) > gap_V_min
                        valid_vel = true;
                        T_Velocity = [V1, V2, V3]; 
                    end
                    vel_retry = vel_retry + 1;
                end
                
                if valid_vel
                    % === C. 计算幅值并验证 ===
                    raw_trend = 1 ./ (T_Range.^2);
                    current_power = mean(raw_trend.^2);
                    true_h_abs = raw_trend / sqrt(current_power);
                    
                    % 检查三个目标中的最大值与最小值之比
                    current_dyn_range = max(true_h_abs) / min(true_h_abs);
                    
                    % 最终卡口：必须 <= 5 (必须与 dataset 生成器完全一样)
                    if current_dyn_range <= max_amp_ratio
                        is_valid = true;
                        true_h = true_h_abs .* exp(1j * 2 * pi * rand(1, num_targets));
                    end
                end
            end
            
            if loop_safety_count > 5000
                R1_min = R1_min * 1.05; % 死循环时略微放宽R1下限
                loop_safety_count = 0;
            end
        end

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

        % --- 4. MTASA 检测算法 (修正后) ---
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

        % 2. 贪心匹配 (不再设置 inf 门限)
        % 我们要确保每个真实目标都分配到一个估计值，即使误差很大
        cost_mat = dist_err_mat;
        matched_pairs = [];

        for pair_idx = 1:num_targets % 强制匹配所有目标
            [min_val, min_linear_idx] = min(cost_mat(:));
            [row, col] = ind2sub(size(cost_mat), min_linear_idx);

            matched_pairs = [matched_pairs; row, col];

            % 记录并移除已匹配的行和列
            cost_mat(row, :) = inf;
            cost_mat(:, col) = inf;
        end

        % 3. 统计误差与检测成功数
        mc_range_sq_err = 0;
        mc_vel_sq_err = 0;
        mc_success_detect = 0; % 用于计算检测概率的计数器

        % 定义检测成功的判定门限（仅用于计算检测概率，不影响 MSE）
        detect_thresh_dist = 1.0 * distance_res;

        for i = 1:num_targets
            t_idx = matched_pairs(i, 1);
            e_idx = matched_pairs(i, 2);

            r_err = dist_err_mat(t_idx, e_idx);
            v_err = vel_err_mat(t_idx, e_idx);

            %% 无论误差多大，全部计入 MSE 累加器
            mc_range_sq_err = mc_range_sq_err + r_err^2;
            mc_vel_sq_err   = mc_vel_sq_err + v_err^2;

            % 如果误差在门限内，才计入“检测成功”
            if r_err < detect_thresh_dist
                mc_success_detect = mc_success_detect + 1;
            end
        end

        % 4. 累加到 SNR 总统计量
        % 注意：此时分母固定为 num_targets，因为所有目标都被计入了
        total_valid_targets_SNR = total_valid_targets_SNR + num_targets;
        sum_sq_err_range = sum_sq_err_range + mc_range_sq_err;
        sum_sq_err_vel   = sum_sq_err_vel + mc_vel_sq_err;

        % 专门给检测概率使用的计数器
        if ~exist('total_success_count_SNR', 'var') || mc == 1
            total_success_count_SNR = 0;
        end
        total_success_count_SNR = total_success_count_SNR + mc_success_detect;

        dist(mc)=mc_range_sq_err;

        %% 打印某个蒙特卡罗下的具体数据
        % num_targets_count = length(true_li);
        % % 打印距离索引信息
        % fprintf('\n--- 距离维度 (Range Index) ---\n');
        % fprintf('%-10s | %-12s | %-12s | %-10s | %-10s\n', '目标', '真实值', '估计值', '偏差', '幅值');
        % for i = 1:num_targets_count
        %     % 匹配对应的估计索引（基于之前 matched_pairs 的逻辑，这里假设 i 是匹配后的顺序）
        %     % 如果只是为了查看当前迭代的原始输出：
        %     fprintf('Path [%d]   | %12.5f | %12.5f | %+.5f | %10.5f\n', ...
        %         i, true_li(i), est_li(i), est_li(i) - true_li(i), abs(true_h(i)));
        % end
        % 
        % % 2. 打印速度维度信息
        % fprintf('\n--- 速度维度 (Velocity Index) ---\n');
        % fprintf('%-10s | %-12s | %-12s | %-10s | %-10s\n', '目标', '真实值', '估计值', '偏差', '幅值');
        % for i = 1:num_targets_count
        %     fprintf('Path [%d]   | %12.5f | %12.5f | %+.5f | %10.5f\n', ...
        %         i, true_ki(i), est_ki(i), est_ki(i) - true_ki(i), abs(true_h(i)));
        % end
        % 
        % figure('Name', 'DD Domain Range-Doppler Heatmap', 'Color', 'w');
        % 
        % % 1. 绘制接收信号幅度热力图
        % % X轴: Doppler (1~N), Y轴: Delay (1~M)
        % imagesc(1:N, 1:M, abs(Y));
        % colormap('jet'); % 使用高对比度颜色
        % colorbar;
        % xlabel('Doppler Index (Column Index)');
        % ylabel('Delay Index (Row I ndex)');
        % title(['Received Signal Y (Magnitude) & Ground Truths (Red X) | SNR=' num2str(current_SNR) 'dB']);
        % hold on;
        % 
        % % 2. 循环绘制真实目标位置
        % for t = 1:num_targets
        %     % --- 坐标转换逻辑 ---
        % 
        %     % A. 距离坐标 (Delay)
        %     % 真实值是 0 ~ M-1，MATLAB索引是 1 ~ M
        %     y_pos = true_li(t) + 1;
        % 
        %     % B. 速度坐标 (Doppler) - 关键！
        %     % 物理真值可能是负数 (例如 -3.5)
        %     % 在矩阵存储中，负频率折叠到了后半部分 (N/2 ~ N)
        %     k_val = true_ki(t);
        % 
        %     if k_val < 0
        %         % 如果是负速度，映射到 N + k_val + 1
        %         % 例如 N=32, k=-2 => index = 32 - 2 + 1 = 31
        %         x_pos = N + k_val + 1;
        %     else
        %         % 如果是正速度，直接 + 1
        %         % 例如 k=2 => index = 2 + 1 = 3
        %         x_pos = k_val + 1;
        %     end
        % 
        %     % 3. 绘制红色的 'x'
        %     plot(x_pos, y_pos, 'rx', 'MarkerSize', 12, 'LineWidth', 2.5);
        % 
        %     % (可选) 添加目标编号标签
        %     text(x_pos, y_pos - 1, sprintf(' T%d', t), ...
        %         'Color', 'white', 'FontWeight', 'bold', 'FontSize', 10);
        % end
        % hold off;
        % % 强制让Y轴方向符合矩阵直观 (行号从上到下增加)
        % set(gca, 'YDir', 'reverse');
        % disp(sqrt(dist(mc)));
        

        
        % % % % figure('Name', 'DD Domain Received Signal (3D Bar)', 'Color', 'w');
        % % % % % 绘制三维柱状图
        % % % % h_bar = bar3(abs(Y));
        % % % % colorbar;
        % % % % % 设置坐标轴标签
        % % % % % 注意：在 MATLAB 矩阵绘图中，行索引对应 Y 轴，列索引对应 X 轴
        % % % % xlabel('Doppler Index (n)', 'FontWeight', 'bold'); % X轴对应列，即 N=32
        % % % % ylabel('Delay Index (m)', 'FontWeight', 'bold');   % Y轴对应行，即 M=64
        % % % % zlabel('Magnitude', 'FontWeight', 'bold');
        % % % % % 设置标题，标出 M 和 N 的值
        % % % % title(sprintf('Received Signal in DD Domain (M=%d, N=%d))', M, N));
        % % % % % 调整视角以便更好地观察导频
        % % % % view(-45, 30); % 方位角 -45度, 仰角 30度
        % % % % % 坐标轴紧凑
        % % % % axis tight;
        

    end
    if total_valid_targets_SNR > 0
        MSE_Range(s_idx) = sqrt(sum_sq_err_range / total_valid_targets_SNR);
        MSE_Velocity(s_idx) = sqrt(sum_sq_err_vel / total_valid_targets_SNR);
    else
        MSE_Range(s_idx) = NaN;
        MSE_Velocity(s_idx) = NaN;
    end

    % 检测概率：成功检测数 / (总实验次数 * 目标数)
    % 更新检测概率：使用专门的成功计数器
    Detection_Prob(s_idx) = total_success_count_SNR / (num_monte_carlo * num_targets);

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
figure('Color', 'w', 'Position', [100, 100, 1200, 400]); % 稍微拉宽窗口

% 1. 距离 MSE 子图
subplot(1,3,1);
semilogy(SNR_vec, MSE_Range, '-ro', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)'); ylabel('RMSE of Range (m)'); % 去掉平方，改为 m
title('Range Estimation Accuracy (RMSE)');


% 2. 速度 MSE 子图
subplot(1,3,2);
semilogy(SNR_vec, MSE_Velocity, '-bs', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)'); ylabel('RMSE of Velocity (m/s)'); % 去掉平方，改为 m/s
title('Velocity Estimation Accuracy (RMSE)');

% 3. 新增：检测概率子图 (非常重要！)
subplot(1,3,3);
plot(SNR_vec, Detection_Prob * 100, '-kP', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)'); ylabel('Detection Probability (%)');
title('Target Detection Performance');
ylim([0 105]); % 纵轴固定在 0-100%




%% ================= 核心算法函数 =================

function [est_li, est_ki] = run_mtasa_detection(Y, dd, M, N, num_targets, delta_f)
% MTASA: Multi-Target Active Sensing Algorithm
% 包含 SIC (串行干扰消除) 初始化 和 PIC (并行干扰消除) 迭代

% 初始化
est_li = zeros(1, num_targets);
est_ki = zeros(1, num_targets);
est_h  = zeros(1, num_targets);

% --- [SIC] 串行干扰消除初始化 ---
% 目的：提供一个较为准确的初值，供 PIC 迭代使用
Residual = Y;
for k = 1:num_targets
    % 1. 两阶段搜索 (包含粗搜+黄金分割精搜)
    [h_hat, l_hat, k_hat] = two_stage_search_matched_golden(Residual, dd, M, N);

    est_li(k) = l_hat;
    est_ki(k) = k_hat;
    est_h(k)  = h_hat;

    % 2. 消除当前检测出的目标
    Signal_Est = generate_echo_physical_truth(dd, h_hat, l_hat, k_hat, M, N);
    Residual = Residual - Signal_Est;
end

% --- [PIC] 并行干扰消除迭代 ---
MaxIter = 10;

for iter = 1:MaxIter
    prev_li = est_li;
    prev_ki = est_ki;

    for k = 1:num_targets
        % 1. 重构干扰 (利用上一轮所有其他目标的参数)
        Interference = zeros(M, N);
        for j = 1:num_targets
            if j ~= k
                Interference = Interference + generate_echo_physical_truth(dd, est_h(j), est_li(j), est_ki(j), M, N);
            end
        end

        % 2. 净化：得到只包含第 k 个目标的信号 r_k
        Y_clean_k = Y - Interference;

        % 3. 【核心修正】严格执行算法 1 (两阶段搜索)
        % 必须重新进行粗搜和精搜，以跳出局部最优
        [h_new, l_new, k_new] = two_stage_search_matched_golden(Y_clean_k, dd, M, N);

        % 4. 更新参数
        est_li(k) = l_new;
        est_ki(k) = k_new;
        est_h(k)  = h_new;
    end

    % 检查收敛 (同时检查距离和速度变化)
    if norm([est_li - prev_li, est_ki - prev_ki]) < 1e-3
        break;
    end
end

% --- 后处理：负速度解卷与距离正数化 ---
mask_neg = est_ki > (N/2);
est_ki(mask_neg) = est_ki(mask_neg) - N;
est_li = mod(est_li, M);
end

% --- 严格遵循论文算法1：两阶段搜索 (含黄金分割) ---
function [h_hat, l_hat, k_hat] = two_stage_search_matched_golden(Y_received, dd, M, N)
% 阶段 1: 粗搜索 (Coarse Search) [cite: 156-163]
% 遍历 DD 域离散网格寻找最大峰值
max_corr = -1; l_coarse_idx = 0; k_coarse_idx = 0;

% 利用 circshift 加速网格遍历,进行lk的粗估计
for l = 0:M-1
    for k = 0:N-1
        dd_shifted = circshift(dd, [l, k]);
        % 论文公式 (18)/(19)
        corr_val = abs(sum(sum(Y_received .* conj(dd_shifted))))^2;
        if corr_val > max_corr
            max_corr = corr_val;
            l_coarse_idx = l; k_coarse_idx = k;
        end
    end
end

% 处理粗搜的负频率索引 (0~N-1 -> -N/2 ~ N/2-1)
if k_coarse_idx > N/2
    k_coarse_idx = k_coarse_idx - N;
end

% 阶段 2: 精搜索 - 二维黄金分割法 (2D Golden Section Search)

% 参数准备
delta_f = 120e3;
T = 1 / delta_f;

% 将粗搜结果转换为物理值 (Delay/Doppler)
tau_c = l_coarse_idx / (M * delta_f);
nu_c  = k_coarse_idx / (N * T);

% 初始搜索区间: [粗搜点 +/- 1个分辨单元]
res_tau = 1 / (M * delta_f);
res_nu  = 1 / (N * T);

a_l = tau_c - res_tau; a_u = tau_c + res_tau; % Delay 区间
b_l = nu_c - res_nu;   b_u = nu_c + res_nu;   % Doppler 区间

mu = (sqrt(5) - 1) / 2; % 黄金分割比例 0.618
Iter = 20; % 黄金分割迭代次数

for i = 1:Iter
    I_a = a_u - a_l;
    I_b = b_u - b_l;

    % 4个试探点
    a1 = a_l + (1 - mu) * I_a; a2 = a_l + mu * I_a;
    b1 = b_l + (1 - mu) * I_b; b2 = b_l + mu * I_b;

    % 计算4个点的相关性，直接根据给定的时延和多普勒，写出对应的回波，然后直接使用回波域接收信号做内积
    vals = zeros(2, 2);
    % (1,1): a1, b1
    E11 = generate_echo_physical_truth_continuous(dd, 1, a1, b1, M, N);
    vals(1,1) = abs(sum(sum(Y_received .* conj(E11))))^2; %如果是计算整体的匹配滤波不能使用这个方法，如果只是算一个点的匹配滤波结果，可以直接做内积
    % (1,2): a1, b2
    E12 = generate_echo_physical_truth_continuous(dd, 1, a1, b2, M, N);
    vals(1,2) = abs(sum(sum(Y_received .* conj(E12))))^2;
    % (2,1): a2, b1
    E21 = generate_echo_physical_truth_continuous(dd, 1, a2, b1, M, N);
    vals(2,1) = abs(sum(sum(Y_received .* conj(E21))))^2;
    % (2,2): a2, b2
    E22 = generate_echo_physical_truth_continuous(dd, 1, a2, b2, M, N);
    vals(2,2) = abs(sum(sum(Y_received .* conj(E22))))^2;

    % 找最大值并缩减区间
    [~, max_idx] = max(vals(:));
    [r_idx, c_idx] = ind2sub([2, 2], max_idx);

    if r_idx == 1 && c_idx == 1      % Max at (a1, b1)
        a_u = a2; b_u = b2;
    elseif r_idx == 1 && c_idx == 2  % Max at (a1, b2)
        a_u = a2; b_l = b1;
    elseif r_idx == 2 && c_idx == 1  % Max at (a2, b1)
        a_l = a1; b_u = b2;
    elseif r_idx == 2 && c_idx == 2  % Max at (a2, b2)
        a_l = a1; b_l = b1;
    end
end

% 最终估计值 (区间中心)
tau_hat_final = (a_l + a_u) / 2;
nu_hat_final  = (b_l + b_u) / 2;

% 转换为连续索引输出
l_hat = tau_hat_final * (M * delta_f);
k_hat = nu_hat_final * (N * T);

% 估计信道系数 h
E_final = generate_echo_physical_truth_continuous(dd, 1, tau_hat_final, nu_hat_final, M, N);
h_hat = sum(sum(Y_received .* conj(E_final))) / sum(sum(abs(E_final).^2));
end

function metric = calculate_metric_matched(Y, dd, l, k, M, N)
E = generate_echo_physical_truth(dd, 1, l, k, M, N);
corr = sum(sum(Y .* conj(E)));
energy = sum(sum(abs(E).^2));
metric = (abs(corr)^2) / energy;
end

% --- 连续值输入 wrapper ---
function Y_out = generate_echo_physical_truth_continuous(dd, h, tau_phy, nu_phy, M, N)
delta_f = 120e3;
T = 1/delta_f;
% 将物理值转换为实数索引 (Float Index)
tau_idx = tau_phy * (M * delta_f);
nu_idx  = nu_phy * (N * T);
Y_out = generate_echo_physical_truth(dd, h, tau_idx, nu_idx, M, N);
end

% --- 核心物理回波生成函数 (支持浮点索引) ---
function Y_out = generate_echo_physical_truth(dd, h, tau_idx, nu_idx, M, N)
% Dirichlet Kernel 物理真值生成
[L_grid, K_grid] = ndgrid(0:M-1, 0:N-1);

% Doppler Kernel (G_term)
k_val = K_grid - nu_idx;
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
Y_out = ifft2(fft2(hw) .* fft2(dd));
end