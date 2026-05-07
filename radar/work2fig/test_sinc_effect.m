% =========================================================================
% MTASA 算法超分辨失败 3D 现场还原 (针对 MC=185 特殊场景) - 中文图表版
% =========================================================================
clc; clear; close all;

%% ================= 1. 参数设置 =================
M = 64;
N = 32;
fc = 64e9;
delta_f = 120e3;
c = 3e8;
delay_res = 1 / (M * delta_f);
doppler_res = 1 / (N * (1/delta_f));
distance_res = c * delay_res / 2;
v_res = doppler_res * c / 2 / fc;

l_max = 20; k_max = 2;
cen_modu = 2; guard_modu = 4;
bits_per_qam_symbol = log2(cen_modu);
bits_per_guard_symbol = log2(guard_modu);
g = 4;

R_max = (l_max/2) * distance_res;
V_max = N * v_res / 2;

num_targets = 3; 
current_SNR = 10; % 固定为 10dB
max_amp_ratio = 5;

guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
data_rows = (l_max + 2):(M - l_max);
guard_rows = [1:(l_max + 1), (M - l_max + 1):M];

fprintf('======================================================\n');
fprintf('启动 MC=185 专属 3D 现场还原...\n');
fprintf('======================================================\n');

%% ================= 2. 生成特定的随机目标 =================
mc = 185; 
rng(mc); % 锁定随机种子

is_valid = false;
while ~is_valid
    % --- A. 生成基础目标 (强制放在前半段，变成 T1 和 T2) ---
    R1 = 2 * distance_res + rand * (R_max / 2 - 2 * distance_res);
    V1 = (rand * 2 - 1) * V_max * 0.7;

    % --- B. 生成伴随目标 (距离和速度相差 0.5 ~ 0.7 个分辨单元) ---
    sign_R = sign(randn); if sign_R == 0, sign_R = 1; end
    sign_V = sign(randn); if sign_V == 0, sign_V = 1; end
    
    delta_R = (0.5 + 0.2 * rand) * distance_res;
    delta_V = (0.5 + 0.2 * rand) * v_res;
    
    R2 = R1 + sign_R * delta_R;
    V2 = V1 + sign_V * delta_V;
    
    if R2 < 2 * distance_res || R2 > (R_max - distance_res) || abs(V2) > (V_max * 0.9)
        continue;
    end

    % --- C. 生成独立的目标 (强制放在后半段，变成 T3) ---
    min_allowed_R_for_T3 = max(R1, R2) + 1.5 * distance_res; 
    max_allowed_R_for_T3 = min(R1, R2) * sqrt(max_amp_ratio); 

    R3_min_bound = max([R_max / 2, min_allowed_R_for_T3]);
    R3_max_bound = min([R_max - distance_res, max_allowed_R_for_T3]);

    if R3_min_bound >= R3_max_bound
        continue; 
    end
    
    R3 = R3_min_bound + rand * (R3_max_bound - R3_min_bound);
    V3 = (rand * 2 - 1) * V_max * 0.9;
    
    if abs(V3 - V1) < 1.0 * v_res || abs(V3 - V2) < 1.0 * v_res
        continue;
    end

    % 整合并排序
    T_Range = [R1, R2, R3];
    T_Velocity = [V1, V2, V3];
    [T_Range, sort_idx] = sort(T_Range);
    T_Velocity = T_Velocity(sort_idx);

    % --- D. 计算幅值并做最终验证 ---
    raw_trend = 1 ./ (T_Range.^2);
    current_power = mean(raw_trend.^2);
    true_h_abs = raw_trend / sqrt(current_power); 

    min_h = min(true_h_abs);
    max_h = max(true_h_abs);
    amp_ratio = max_h / min_h;

    if min_h > 0.0065 && amp_ratio < max_amp_ratio
        is_valid = true;
        true_h = true_h_abs .* exp(1j * 2 * pi * rand(1, num_targets));
    end
end

true_li = T_Range ./ distance_res;
true_ki = T_Velocity ./ v_res;

%% ================= 3. 生成发射与接收信号 =================
dd = zeros(M, N);
dd(1,1) = 1000; 

num_groups_data = N / g;
for row = data_rows
    for j = 1:num_groups_data
        im_bits = randi([0 1], 1, log2(nchoosek(g, 1)));
        qam_bits = randi([0 1], 1, bits_per_qam_symbol);
        active_idx = bi2de(im_bits, 'left-msb') + 1;
        sym_val = qammod(bi2de(qam_bits, 'left-msb'), cen_modu, 'UnitAveragePower', true);
        dd(row, (j-1)*g + active_idx) = sym_val; 
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
            sym_val = qammod(bi2de(randi([0 1], 1, bits_per_guard_symbol), 'left-msb'), guard_modu, 'UnitAveragePower', true);
            dd(row, guard_cols_range(i_col + act_idx - 1)) = sym_val; 
        end
        i_col = i_col + grp_sz;
    end
end

current_signal_power = sum(abs(dd(:)).^2) / (M * N);
target_power = 10^(current_SNR / 10);
dd = dd * sqrt(target_power / current_signal_power);

Y_clean = zeros(M, N);
for p = 1:num_targets
    Y_clean = Y_clean + generate_echo_physical_truth(dd, true_h(p), true_li(p), true_ki(p), M, N);
end
noise = sqrt(1/2) * (randn(M, N) + 1j * randn(M, N));
Y = Y_clean + noise;

%% ================= 4. 运行 MTASA (带 3D 监控) =================
[est_li, est_ki, est_h] = run_mtasa_detection_with_monitor(Y, dd, M, N, num_targets, delta_f, true_li, true_ki);

%% ================= 5. 计算偏差并打印 =================
fprintf('\n--- 目标匹配与估计结果 (MC=185 专属还原) ---\n');
fprintf('%-11s | %-27s | %-27s | %-27s\n', '配对状态', '距离维度 (Range Index)', '速度维度 (Velocity Index)', '增益幅度 (Gain Amplitude)');
fprintf('%-11s | %-8s %-8s %-8s | %-8s %-8s %-8s | %-8s %-8s %-8s\n', ...
    'T <-> E', '真实值', '估计值', '偏差', '真实值', '估计值', '偏差', '真实值', '估计值', '偏差');

for t_idx = 1:num_targets
    e_idx = t_idx;
    r_true = true_li(t_idx); r_est = est_li(e_idx); r_diff = r_est - r_true;
    v_true = true_ki(t_idx); v_est = est_ki(e_idx); v_diff = v_est - v_true;
    
    if abs(v_diff) > N/2
        v_diff = v_diff - sign(v_diff)*N;
    end
    
    h_true_amp = abs(true_h(t_idx)); h_est_amp = abs(est_h(e_idx)); h_diff = h_est_amp - h_true_amp;
    
    if abs(r_diff) > 1.0
        pair_str = sprintf('T[%d]<!>E[%d]', t_idx, e_idx);
    else
        pair_str = sprintf('T[%d]<->E[%d]', t_idx, e_idx);
    end
    
    fprintf('%-11s | %8.5f %8.5f %+8.5f | %8.5f %8.5f %+8.5f | %8.5f %8.5f %+8.5f\n', ...
        pair_str, r_true, r_est, r_diff, v_true, v_est, v_diff, h_true_amp, h_est_amp, h_diff);
end


%% ======================= 核心函数区 =======================

function [est_li, est_ki, est_h] = run_mtasa_detection_with_monitor(Y, dd, M, N, num_targets, delta_f, true_li, true_ki)
    est_li = zeros(1, num_targets);
    est_ki = zeros(1, num_targets);
    est_h  = zeros(1, num_targets);

    % --- [SIC] 初始化 ---
    Residual = Y;
    for k = 1:num_targets
        
        % ====== 【3D 监控摄像头：捕捉寻找 T3 时的残余能量】 ======
        if k == 3 
            fprintf('\n>> [监控] 正在生成寻找 T3 时的 2D 残差能量曲面图...\n');
            Corr_Map = zeros(M, N);
            for l_idx = 0:M-1
                for k_idx = 0:N-1
                    dd_shifted = circshift(dd, [l_idx, k_idx]);
                    Corr_Map(l_idx+1, k_idx+1) = abs(sum(sum(Residual .* conj(dd_shifted))));
                end
            end
            
            figure('Name', '2D 残差能量曲面图 (搜索 T3)', 'Position', [100, 100, 900, 600], 'Color', 'w');
            surf(0:N-1, 0:M-1, Corr_Map, 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'DisplayName', '残差能量底图');
            colormap('jet'); colorbar; hold on;
            
            [max_val, max_idx] = max(Corr_Map(:));
            [max_l, max_k] = ind2sub([M, N], max_idx);
            
            t3_k_disp = mod(true_ki(3), N); 
            t3_l_disp = true_li(3);
            t3_val = Corr_Map(round(t3_l_disp)+1, round(t3_k_disp)+1);
            
            % 标记算法眼中的“最高峰”(幻影)
            plot3(max_k-1, max_l-1, max_val, 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'r', 'DisplayName', '幻影峰 (算法将被诱导抓取此处)');
            
            % 标记真实的 T3 躲在哪儿
            plot3(t3_k_disp, t3_l_disp, t3_val, 'gp', 'MarkerSize', 18, 'MarkerFaceColor', 'g', 'DisplayName', '真实 T3 位置 (被垃圾能量掩没)');
            
            % 设置中文字体，防止乱码 (优先使用微软雅黑或黑体)
            title('剔除 T1 和 T2 后的 2D 残差能量景观 (准备搜索 T3)', 'FontSize', 14, 'FontName', 'Microsoft YaHei');
            xlabel('多普勒索引 (速度)', 'FontSize', 12, 'FontName', 'Microsoft YaHei'); 
            ylabel('时延索引 (距离)', 'FontSize', 12, 'FontName', 'Microsoft YaHei'); 
            zlabel('相关性能量幅值', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
            legend('Location', 'northeast', 'FontSize', 12, 'FontName', 'Microsoft YaHei');
            
            % 设置坐标轴的字体
            set(gca, 'FontName', 'Microsoft YaHei', 'FontSize', 10);
            
            view(-40, 30); 
            grid on; hold off;
            fprintf('>> [监控完成] 图表已生成。继续执行算法...\n\n');
        end
        % ========================================================

        [h_hat, l_hat, k_hat] = two_stage_search_matched_golden(Residual, dd, M, N);
        est_li(k) = l_hat;
        est_ki(k) = k_hat;
        est_h(k)  = h_hat;
        Signal_Est = generate_echo_physical_truth(dd, h_hat, l_hat, k_hat, M, N);
        Residual = Residual - Signal_Est;
    end

    % --- [PIC] 并行干扰消除迭代 ---
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

    delta_f = 120e3; T = 1 / delta_f;
    tau_c = l_coarse_idx / (M * delta_f);
    nu_c  = k_coarse_idx / (N * T);

    res_tau = 1 / (M * delta_f); res_nu  = 1 / (N * T);

    a_l = tau_c - 1.0 * res_tau; a_u = tau_c + 1.0 * res_tau;
    b_l = nu_c - 1.0 * res_nu;   b_u = nu_c + 1.0 * res_nu;

    mu = (sqrt(5) - 1) / 2; Iter = 20;
    for i = 1:Iter
        I_a = a_u - a_l; I_b = b_u - b_l;
        a1 = a_l + (1 - mu) * I_a; a2 = a_l + mu * I_a;
        b1 = b_l + (1 - mu) * I_b; b2 = b_l + mu * I_b;

        vals = zeros(2, 2);
        vals(1,1) = abs(sum(sum(Y_received .* conj(generate_echo_physical_truth_continuous(dd, 1, a1, b1, M, N)))))^2;
        vals(1,2) = abs(sum(sum(Y_received .* conj(generate_echo_physical_truth_continuous(dd, 1, a1, b2, M, N)))))^2;
        vals(2,1) = abs(sum(sum(Y_received .* conj(generate_echo_physical_truth_continuous(dd, 1, a2, b1, M, N)))))^2;
        vals(2,2) = abs(sum(sum(Y_received .* conj(generate_echo_physical_truth_continuous(dd, 1, a2, b2, M, N)))))^2;

        [~, max_idx] = max(vals(:));
        [r_idx, c_idx] = ind2sub([2, 2], max_idx);
        if r_idx == 1 && c_idx == 1, a_u = a2; b_u = b2;
        elseif r_idx == 1 && c_idx == 2, a_u = a2; b_l = b1;
        elseif r_idx == 2 && c_idx == 1, a_l = a1; b_u = b2;
        elseif r_idx == 2 && c_idx == 2, a_l = a1; b_l = b1;
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
    delta_f = 120e3; T = 1/delta_f;
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