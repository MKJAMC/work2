%% ================= OTFS 普适多目标训练集生成 (物理衰减 1/R^2 升级版) =================
% 核心特性:
% 1. 目标数随机 (1~5)
% 2. SNR 连续随机 (10~30dB)
% 3. 三态混合拓扑 (完全分散 / 完全扎堆 / 随机部分扎堆)
% 4. 严格遵守物理学大尺度衰落: 幅值与距离的平方成反比 (1/R^2)
% 5. 空间域预限制: 保证自然衰减后的增益动态范围 < 5 倍

clc; clear; close all;

% ================= 1. 基础物理参数 =================
M = 64; N = 32; fc = 64e9; delta_f = 120e3; c = 3e8;
distance_res = c / (2 * M * delta_f);
v_res = (c * delta_f) / (2 * N * fc);

num_samples = 8000; % 生成样本数
P_max = 5;          % 最大支持 5 个目标

l_max = 20; k_max = 2;
R_max = (l_max/2) * distance_res;
V_max = (N/2) * v_res;

% --- IM-OTFS 调制参数 ---
g = 4; 
cen_modu = 2; guard_modu = 4; 
bits_per_qam = log2(cen_modu);
bits_per_guard = log2(guard_modu);

guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
data_rows = (l_max + 2):(M - l_max);
guard_rows = [1:(l_max + 1), (M - l_max + 1):M];

% --- 预分配存储 ---
Y_data = zeros(M, N, num_samples, 'like', 1j);
X_data = zeros(M, N, num_samples, 'like', 1j); 
Labels = zeros(num_samples, P_max, 4); 
SNR_Index = zeros(num_samples, 1);

fprintf('开始生成 普适版 多目标训练数据 (基于 1/R^2 物理距离衰减)...\n');
timer_start = tic;

for mc = 1:num_samples
    if mod(mc, 500) == 0
        fprintf('已生成 %d / %d 样本...\n', mc, num_samples); 
    end
    
    current_SNR = 10 + rand() * 20; 
    SNR_Index(mc) = current_SNR;
    
    power_persym = 10.^(current_SNR/10);
    ans2 = (M - 2*l_max - 1)*N + (N - 4*k_max - 1)*(2*l_max + 1);
    power = power_persym * ans2;
    factor0 = sqrt(power ./ 430);
    
    % =========================================================
    % B. 随机目标位置生成 (引入 L_max / L_min < 2.2 物理约束)
    % =========================================================
    scenario_type = rand();
    
    % 为了保证 (远/近)^2 < 5，时延网格的比例不能超过 sqrt(5) 约等于 2.23
    % 设定一个基准锚点，所有生成的点必须在这个比例限制内
    base_li = 2.0 + rand() * 4.0; % 锚点在 2.0 到 6.0 之间
    max_li_allow = min(base_li * 2.2, l_max); 
    
    if scenario_type < 0.33
        % [场景 1] 均匀分散
        P = randi([1, P_max]); 
        true_li = zeros(1, P); true_ki = zeros(1, P);
        for p = 1:P
            is_valid = false;
            while ~is_valid
                cand_li = base_li + rand() * (max_li_allow - base_li);
                cand_ki = (rand() * 2 - 1) * (V_max * 0.8)/v_res;
                if p == 1, is_valid = true;
                else
                    dist = sqrt((true_li(1:p-1) - cand_li).^2 + (true_ki(1:p-1) - cand_ki).^2);
                    if min(dist) > 2.0, is_valid = true; end
                end
            end
            true_li(p) = cand_li; true_ki(p) = cand_ki;
        end
        
    elseif scenario_type < 0.66
        % [场景 2] 紧凑扎堆 (0.3 ~ 0.6 bins)
        P = randi([2, 3]); 
        true_li = zeros(1, P); true_ki = zeros(1, P);
        
        center_li = base_li + rand() * (max_li_allow - base_li - 1.0);
        center_ki = (rand() * 2 - 1) * (V_max * 0.5)/v_res;
        true_li(1) = center_li; true_ki(1) = center_ki;
        
        for p = 2:P
            is_valid = false; attempts = 0;
            while ~is_valid && attempts < 200
                sign_l = randi([0 1])*2 - 1; sign_k = randi([0 1])*2 - 1; 
                cand_li = center_li + sign_l * (0.3 + rand() * 0.3); 
                cand_ki = center_ki + sign_k * (0.3 + rand() * 0.3); 
                
                dist_l = abs(true_li(1:p-1) - cand_li);
                dist_k = abs(true_ki(1:p-1) - cand_ki);
                dist_euclid = sqrt(dist_l.^2 + dist_k.^2);
                
                if min(dist_l) >= 0.3 && min(dist_k) >= 0.3 && max(dist_euclid) <= 1.2
                    % 额外约束：必须在物理允许的动态框内
                    if cand_li >= base_li && cand_li <= max_li_allow
                        is_valid = true; 
                    end
                end
                attempts = attempts + 1;
            end
            if ~is_valid % 兜底
                cand_li = center_li + 0.4 * (randi([0 1]*2-1));
                cand_ki = center_ki + 0.4 * (randi([0 1]*2-1));
            end
            true_li(p) = cand_li; true_ki(p) = cand_ki;
        end
        
    else
        % [场景 3] 混合状态 (扎堆簇 + 远处散点)
        cluster_size = randi([2, 3]); scatter_size = randi([1, 2]);   
        P = cluster_size + scatter_size; 
        true_li = zeros(1, P); true_ki = zeros(1, P);
        
        center_li = base_li + rand() * ((max_li_allow - base_li) * 0.4);
        center_ki = (rand() * 2 - 1) * (V_max * 0.5)/v_res;
        true_li(1) = center_li; true_ki(1) = center_ki;
        
        for p = 2:P
            is_valid = false; attempts = 0;
            while ~is_valid && attempts < 200
                if p <= cluster_size
                    % 簇内目标
                    sign_l = randi([0 1])*2 - 1; sign_k = randi([0 1])*2 - 1;
                    cand_li = center_li + sign_l * (0.3 + rand() * 0.3);
                    cand_ki = center_ki + sign_k * (0.3 + rand() * 0.3);
                else
                    % 散点：要在簇外(>2.0)，且在 max_li_allow 之内
                    cand_li = center_li + 2.0 + rand() * (max_li_allow - center_li - 2.0);
                    cand_ki = (rand() * 2 - 1) * (V_max * 0.8)/v_res;
                end
                
                dist_l = abs(true_li(1:p-1) - cand_li);
                dist_k = abs(true_ki(1:p-1) - cand_ki);
                dist_euclid = sqrt(dist_l.^2 + dist_k.^2);
                
                if p <= cluster_size
                    if min(dist_l) >= 0.3 && min(dist_k) >= 0.3 && max(dist_euclid) <= 1.2
                        is_valid = true; 
                    end
                else
                    if min(dist_euclid) >= 2.0 && cand_li >= base_li && cand_li <= max_li_allow
                        is_valid = true; 
                    end
                end
                attempts = attempts + 1;
            end
            if ~is_valid && p > cluster_size % 散点兜底
                cand_li = min(center_li + 2.5, max_li_allow);
                cand_ki = 2.0 + rand() * 5;
            end
            true_li(p) = cand_li; true_ki(p) = cand_ki;
        end
    end
    
    T_Range = true_li * distance_res;
    T_Velocity = true_ki * v_res;

    % =========================================================
    % C. 物理信道幅值生成
    % =========================================================
    % 1. 物理距离衰减 (Radar方程: 振幅与 R^2 成反比)
    distance_attenuation = 1 ./ (T_Range.^2);
    
    % 2. 【修改】直接使用物理衰减，彻底移除瑞利随机快衰落
    % rayleigh_fading = abs(randn(1, P) + 1j * randn(1, P)); 
    raw_h_abs = distance_attenuation; 
    
    % 3. 防深衰落保护 (由于移除了瑞利衰落，这一步其实很少会触发，但保留作为安全兜底)
    raw_h_abs = raw_h_abs / max(raw_h_abs); 
    raw_h_abs = max(raw_h_abs, 0.15); 
    
    % 4. 平均功率归一化 (保证发射端 SNR 绝对精准)
    current_mean_power = mean(raw_h_abs.^2); 
    true_h_abs = raw_h_abs / sqrt(current_mean_power); 
    
    % 5. 附加随机相位
    true_h = true_h_abs .* exp(1j * 2 * pi * rand(1, P));
    
    % =========================================================
    % D & E. 信号合成与加噪保持不变
    % =========================================================
    dd = zeros(M, N); dd(1,1) = 1000; % Pilot
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
    
    Y_clean = zeros(M, N);
    for p = 1:P
        Y_clean = Y_clean + generate_echo_physical_truth(dd, true_h(p), true_li(p), true_ki(p), M, N);
    end
    noise = (randn(M, N) + 1j * randn(M, N)) * sqrt(0.5);
    
    Y_data(:,:,mc) = Y_clean + noise;
    X_data(:,:,mc) = dd; 
    
    for p = 1:P
        Labels(mc, p, :) = [true_li(p), true_ki(p), abs(true_h(p)), angle(true_h(p))];
    end
end
toc(timer_start);

% ================= 4. 导出为 H5 =================
% 指定你的科研数据保存路径
save_dir = 'D:\研究生\科研\frac_doppler_channel_est\radar\radar_dection\data_h5';

% 如果文件夹不存在，则自动创建，防止最后一步报错白跑
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('已自动创建文件夹: %s\n', save_dir);
end

% 拼接完整的文件路径
filename = fullfile(save_dir, 'OTFS_Train_Set_Physical.h5'); 

if exist(filename, 'file')
    delete(filename); 
end

h5create(filename, '/Y_real', size(real(Y_data))); h5write(filename, '/Y_real', real(Y_data));
h5create(filename, '/Y_imag', size(imag(Y_data))); h5write(filename, '/Y_imag', imag(Y_data));
h5create(filename, '/X_real', size(real(X_data))); h5write(filename, '/X_real', real(X_data));
h5create(filename, '/X_imag', size(imag(X_data))); h5write(filename, '/X_imag', imag(X_data));
h5create(filename, '/labels', size(Labels)); h5write(filename, '/labels', Labels);
h5create(filename, '/snr_index', size(SNR_Index)); h5write(filename, '/snr_index', SNR_Index);

fprintf('训练集生成完毕！共 %d 个样本，已保存至: %s\n', num_samples, filename);

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