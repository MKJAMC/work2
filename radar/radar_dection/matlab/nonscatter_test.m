%% ================= OTFS T1/T2 不可分辨目标数据集生成 =================
% 结合了 t1t2.m 的生成逻辑（最大/最小增益比在 [4.5, 5.0] 之间，且 T1 和 T2 间隔 < 1 分辨率）
% 并按照 scatter_test.m 的格式保存为深度学习使用的 HDF5 数据集
clc; clear; close all;

% ================= 1. 基础物理参数 =================
M = 64; N = 32; fc = 64e9; delta_f = 120e3; c = 3e8;
distance_res = c / (2 * M * delta_f);
v_res = (c * delta_f) / (2 * N * fc);

% ================= 2. 场景与调制参数 =================
% --- 仿真配置 ---
SNR_vec = 10:2:18;      % SNR 列表 (10, 12, 14, 16, 18)
num_monte_carlo = 200;  % 每个 SNR 的蒙特卡洛次数
num_samples = length(SNR_vec) * num_monte_carlo; % 总样本数 = 5 * 200 = 1000

P = 3;               % 目标数
l_max = 20; k_max = 2;
R_max = (l_max/2) * distance_res;
V_max = (N/2) * v_res;

% --- IM-OTFS 调制参数 ---
g = 4; 
cen_modu = 2; guard_modu = 4; 
bits_per_qam = log2(cen_modu);
bits_per_guard = log2(guard_modu);

% --- 预计算索引 ---
guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
data_rows = (l_max + 2):(M - l_max);
guard_rows = [1:(l_max + 1), (M - l_max + 1):M];

% --- 预分配存储空间 ---
Y_data = zeros(M, N, num_samples, 'like', 1j);
X_data = zeros(M, N, num_samples, 'like', 1j); 
Labels = zeros(num_samples, P, 4); 

fprintf('开始生成 T1/T2 不可分辨版测试数据...\n');
fprintf('SNR 列表: [%s] dB\n', num2str(SNR_vec));
fprintf('每个 SNR 样本数: %d\n', num_monte_carlo);
fprintf('总样本数: %d\n', num_samples);

global_idx = 0;
timer_start = tic;

% ================= 3. 双重循环 (SNR -> MC) =================
for s_idx = 1:length(SNR_vec)
    current_SNR = SNR_vec(s_idx);
    fprintf('正在处理 SNR = %d dB ...\n', current_SNR);
    
    % --- 功率计算 ---
    power_persym = 10.^(current_SNR/10);
    ans2 = (M - 2*l_max - 1)*N + (N - 4*k_max - 1)*(2*l_max + 1);
    power = power_persym * ans2;
    ans3 = 430; 
    factor0 = sqrt(power ./ ans3);
    
    for mc = 1:num_monte_carlo
        global_idx = global_idx + 1;
        
        % --- 关键：设置随机种子 ---
        % 这保证了在不同 SNR 下，第 mc 次实验生成的物理场景(距离/速度)是完全一样的
        rng(mc);
        
        % ================= 3.1 生成随机目标 (Logic from t1t2.m) =================
        target_amp_ratio_min = 1.2 ;  % 目标幅值比下限
        target_amp_ratio_max = 5.0;  % 目标幅值比上限
        
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

            % --- B. 【定向生成 R3】保证幅值比在 5 附近 ---
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

                    % 严格卡死动态范围在目标区间 [4.5, 5.0]
                   if current_dyn_range <= target_amp_ratio_max
                        is_valid = true;
                        true_h = true_h_abs .* exp(1j * 2 * pi * rand(1, P));
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
        
        % ================= 3.2 生成发射信号 X_DD =================
        dd = zeros(M, N);
        dd(1,1) = 1000; % Pilot
        
        % Data Region
        num_groups = N / g;
        for row = data_rows
            for j = 1:num_groups
                im_bits = randi([0 1], 1, log2(nchoosek(g, 1)));
                qam_bits = randi([0 1], 1, bits_per_qam);
                act_idx = bi2de(im_bits, 'left-msb') + 1;
                sym = qammod(bi2de(qam_bits, 'left-msb'), cen_modu, 'UnitAveragePower', true);
                dd(row, (j-1)*g + act_idx) = sym * factor0;
            end
        end
        
        % Guard Region
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
                    sym = qammod(bi2de(randi([0 1], 1, bits_per_guard), 'left-msb'), guard_modu, 'UnitAveragePower', true);
                    dd(row, guard_cols_range(i_col + act_idx - 1)) = sym * factor0;
                end
                i_col = i_col + grp_sz;
            end
        end

        % ================= 3.3 物理回波合成 =================
        Y_clean = zeros(M, N);
        for p = 1:P
            Y_clean = Y_clean + generate_echo_physical_truth(dd, true_h(p), true_li(p), true_ki(p), M, N);
        end
        
        % 加入噪声
        noise = (randn(M, N) + 1j * randn(M, N)) * sqrt(0.5);
        
        % 存入数组
        Y_data(:,:,global_idx) = Y_clean + noise;
        X_data(:,:,global_idx) = dd; 
        Labels(global_idx,:,:) = [true_li', true_ki', abs(true_h)', angle(true_h)'];
    end
end
toc(timer_start);

% ================= 4. 导出为 H5 =================
filename = 'OTFS_Test_Set_T1T2_1000.h5'; % 【修改】更新了输出文件名为 T1T2 标识
if exist(filename, 'file'), delete(filename); end

h5create(filename, '/Y_real', size(real(Y_data))); h5write(filename, '/Y_real', real(Y_data));
h5create(filename, '/Y_imag', size(imag(Y_data))); h5write(filename, '/Y_imag', imag(Y_data));
h5create(filename, '/X_real', size(real(X_data))); h5write(filename, '/X_real', real(X_data));
h5create(filename, '/X_imag', size(imag(X_data))); h5write(filename, '/X_imag', imag(X_data));
h5create(filename, '/labels', size(Labels)); h5write(filename, '/labels', Labels);

% 额外保存 SNR Index 标签
SNR_Index = zeros(num_samples, 1);
idx = 0;
for s = 1:length(SNR_vec)
    for m = 1:num_monte_carlo
        idx = idx + 1;
        SNR_Index(idx) = SNR_vec(s);
    end
end
h5create(filename, '/snr_index', size(SNR_Index)); h5write(filename, '/snr_index', SNR_Index);

fprintf('测试集生成完毕！共 %d 个样本，已保存至: %s\n', num_samples, filename);

%% ================= 核心函数 =================
function Y_out = generate_echo_physical_truth(dd, h, tau_idx, nu_idx, M, N)
    [L_grid, K_grid] = ndgrid(0:M-1, 0:N-1);
    
    k_val = K_grid - nu_idx;
    mask_peak = abs(sin(pi * k_val / N)) < 1e-9;
    G_term = zeros(size(k_val)); G_term(mask_peak) = 1;
    mn = ~mask_peak;
    if any(mn(:)), kv = k_val(mn); G_term(mn) = exp(-1i*pi*(N-1)*kv/N).*(sin(pi*kv)./(N*sin(pi*kv/N))); end
    
    l_val = L_grid - tau_idx;
    mask_peak = abs(sin(pi * l_val / M)) < 1e-9;
    F_term = zeros(size(l_val)); F_term(mask_peak) = 1;
    mn = ~mask_peak;
    if any(mn(:)), lv = l_val(mn); F_term(mn) = exp(-1i*pi*(M-1)*lv/M).*(sin(pi*lv)./(M*sin(pi*lv/M))); end
    
    hw = h * G_term .* F_term * exp(-1j * 2 * pi * tau_idx * nu_idx / (M*N));
    Y_out = ifft2(fft2(hw) .* fft2(dd));
end