%% ================= OTFS 分散目标数据集生成 (基于 Scatter 逻辑) =================
clc; clear; close all;

% ================= 1. 基础物理参数 =================
M = 64; N = 32; fc = 64e9; delta_f = 120e3; c = 3e8;
distance_res = c / (2 * M * delta_f);
v_res = (c * delta_f) / (2 * N * fc);

% ================= 2. 场景与调制参数 =================
% --- 仿真配置 ---
SNR_vec = 10:2:18;      % 【修改点】SNR 列表 (10, 12, 14, 16, 18)
num_monte_carlo = 200;  % 【修改点】每个 SNR 的蒙特卡洛次数
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

fprintf('开始生成严格分散版测试数据 (Scatter Logic)...\n');
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
        
        % ================= 3.1 生成随机目标 (Logic from scatter.m) =================
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
                    
                    % 最终卡口：必须 < 5
                    if current_dyn_range <= max_amp_ratio
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
filename = 'OTFS_Test_Set_Scatter_1000.h5'; % 【修改】更新了输出文件名
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