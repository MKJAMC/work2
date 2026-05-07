% 通过调整true_h最大与最小之间的差距不超过x db
% PIC的迭代次数 10
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
num_targets = 3; % 目标个数
SNR_vec = 10:2:18;  % SNR 范围 (dB)
num_monte_carlo = 10; % 蒙特卡洛次数 (测试时建议先改小，比如10，跑通后再改回500)
max_amp_ratio = 5;

% 预计算索引参数
guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
data_rows = (l_max + 2):(M - l_max);
guard_rows = [1:(l_max + 1), (M - l_max + 1):M];

% 结果存储
MSE_Range = zeros(length(SNR_vec), 1);
MSE_Velocity = zeros(length(SNR_vec), 1);
Detection_Prob = zeros(length(SNR_vec), 1); % 新增初始化

fprintf('======================================================\n');
fprintf('开始 Monte Carlo 仿真 (终极修正版：动态功率 + 严格配对)\n');
fprintf('SNR 范围: [%s] dB\n', num2str(SNR_vec));
fprintf('======================================================\n');

%% ================= 主循环 (SNR) =================
total_timer = tic;

for s_idx = 1:length(SNR_vec)
    current_SNR = SNR_vec(s_idx);
    fprintf('\n正在仿真 SNR = %d dB ...\n', current_SNR);

    % 累加误差变量
    total_valid_targets_SNR = 0; % 关键：每个 SNR 点开始前，有效目标计数清零
    sum_sq_err_range = 0;
    sum_sq_err_vel = 0;
    dist = zeros(1, num_monte_carlo);

    %% 蒙特卡洛循环
    for mc = 1:num_monte_carlo
        fprintf('MC: %d/%d\n', mc, num_monte_carlo);
        mc=185;
        rng(mc); % 固定随机种子，方便复现

        %% === 约束条件 ===
        % 物理间隔
        gap_R = 0.5 * distance_res;
        gap_V = 0.5 * v_res;

        %% === 1. 生成随机目标 (增加幅值动态范围限制) ===
        is_valid = false;
        loop_safety_count = 0;
        while ~is_valid
            loop_safety_count = loop_safety_count + 1;
            
            % --- A. 生成基础目标 (T1T2离得近) ---
            R1 = 2 * distance_res + rand * (R_max / 2 - 2 * distance_res);
            V1 = (rand * 2 - 1) * V_max * 0.7;

            % --- B. 生成伴随目标 (强制与 R1 的时延和多普勒相差 0.5 ~ 0.7 个分辨单元) ---
            sign_R = sign(randn); if sign_R == 0, sign_R = 1; end
            sign_V = sign(randn); if sign_V == 0, sign_V = 1; end

            % 锁死 0.5 到 0.7 的差距
            delta_R = (0.5 + 0.2 * rand) * distance_res;
            delta_V = (0.5 + 0.2 * rand) * v_res;

            R2 = R1 + sign_R * delta_R;
            V2 = V1 + sign_V * delta_V;

            if R2 < 2 * distance_res || R2 > (R_max - distance_res) || abs(V2) > (V_max * 0.9)
                continue;
            end

            % --- C. 生成独立的目标 (强制放在后半段，排序后变成 T3) ---
            % 由于 T1/T2 现在距离近（信号强），T3 距离远（信号弱）
            % 我们需要严格控制 T3 不能太远，否则 T1/T2 和 T3 的增益比会超过 max_amp_ratio (5)

            % 1. 物理安全距离限制：至少和双胞胎隔开 1.5 个分辨单元
            min_allowed_R_for_T3 = max(R1, R2) + 1.5 * distance_res; 

            % 2. 增益比限制：由于幅度与 1/R^2 成正比，保证最强(近)与最弱(远)的幅度比 < max_amp_ratio
            max_allowed_R_for_T3 = min(R1, R2) * sqrt(max_amp_ratio); 

            % R3 的下界：必须在后半段，且满足最小物理间隔
            R3_min_bound = max([R_max / 2, min_allowed_R_for_T3]);
            % R3 的上界：不能超过系统最大距离，且不能违背增益比限制
            R3_max_bound = min([R_max - distance_res, max_allowed_R_for_T3]);

            if R3_min_bound >= R3_max_bound
                continue; % 无法同时满足边界条件（比如T1/T2太近导致T3按比例不能放太远），重新生成
            end

            R3 = R3_min_bound + rand * (R3_max_bound - R3_min_bound);
            V3 = (rand * 2 - 1) * V_max * 0.9;

            % 确保独立目标在速度上与双胞胎至少间隔 1 个分辨单元
            if abs(V3 - V1) < 1.0 * v_res || abs(V3 - V2) < 1.0 * v_res
                continue;
            end

            % T2T3离得近
            % --- A. 生成基础目标 (强制放在后半段) ---
            % % 删除了 50% 的随机概率，强制 R1 在较远的位置
            % R1 = R_max / 2 + rand * (R_max / 2 - 2 * distance_res);
            % V1 = (rand * 2 - 1) * V_max * 0.7;
            % 
            % % --- B. 生成伴随目标 (强制与 R1 的时延和多普勒相差 0.5 ~ 0.7 个分辨单元) ---
            % sign_R = sign(randn); if sign_R == 0, sign_R = 1; end
            % sign_V = sign(randn); if sign_V == 0, sign_V = 1; end
            % % 锁死 0.5 到 0.7 的差距
            % delta_R = (0.5 + 0.2 * rand) * distance_res;
            % delta_V = (0.5 + 0.2 * rand) * v_res;
            % R2 = R1 + sign_R * delta_R;
            % V2 = V1 + sign_V * delta_V;
            % 
            % if R2 < 2 * distance_res || R2 > (R_max - distance_res) || abs(V2) > (V_max * 0.9)
            %     continue;
            % end
            % 
            % % --- C. 生成独立的目标 (强制放在前半段，变成 T1) ---
            % min_allowed_R_for_T3 = max(R1, R2) / sqrt(max_amp_ratio);
            % max_allowed_R_for_T3 = min(R1, R2) * sqrt(max_amp_ratio);
            % 
            % R3_min_bound = max([2 * distance_res, min_allowed_R_for_T3]);
            % % 【关键修改】强制 R3 最大不能超过 R1 和 R2，且保持至少 1.5 个分辨单元的安全距离
            % R3_max_bound = min([R_max - distance_res, max_allowed_R_for_T3, min(R1, R2) - 1.5 * distance_res]);
            % 
            % if R3_min_bound >= R3_max_bound
            %     continue; % 无法同时满足边界条件，重新生成
            % end
            % 
            % R3 = R3_min_bound + rand * (R3_max_bound - R3_min_bound);
            % V3 = (rand * 2 - 1) * V_max * 0.9;
            % 
            % % 确保独立目标在速度上与双胞胎至少间隔 1 个分辨单元
            % if abs(V3 - V1) < 1.0 * v_res || abs(V3 - V2) < 1.0 * v_res
            %     continue;
            % end
            

            % 整合并排序
            T_Range = [R1, R2, R3];
            T_Velocity = [V1, V2, V3];
            [T_Range, sort_idx] = sort(T_Range);
            T_Velocity = T_Velocity(sort_idx);



            % T1T2T3随机距离
            % is_valid = false;
            % loop_safety_count = 0;
            %
            % while ~is_valid
            %     loop_safety_count = loop_safety_count + 1;
            %
            %     --- A. 生成目标 T1 (核心修改：50%概率前半段，50%概率后半段) ---
            %     if rand > 0.5
            %         放在前半段：T3大概率会生成在右侧。排序后，这对双胞胎会变成 T1 和 T2。
            %         R1 = 2 * distance_res + rand * (R_max / 2 - 2 * distance_res);
            %     else
            %         放在后半段：T3大概率会生成在左侧。排序后，T3变成T1，这对双胞胎会变成 T2 和 T3。
            %         R1 = R_max / 2 + rand * (R_max / 2 - 2 * distance_res);
            %     end
            %
            %     V1 = (rand * 2 - 1) * V_max * 0.7;
            %
            %     --- B. 生成目标 T2 (强制与 T1 的时延和多普勒相差 0.5 ~ 0.7 个分辨单元) ---
            %     sign_R = sign(randn); if sign_R == 0, sign_R = 1; end
            %     sign_V = sign(randn); if sign_V == 0, sign_V = 1; end
            %
            %     锁死 0.5 到 0.7 的差距
            %     delta_R = (0.5 + 0.2 * rand) * distance_res;
            %     delta_V = (0.5 + 0.2 * rand) * v_res;
            %
            %     R2 = R1 + sign_R * delta_R;
            %     V2 = V1 + sign_V * delta_V;
            %
            %     if R2 < 2 * distance_res || R2 > (R_max - distance_res) || abs(V2) > (V_max * 0.9)
            %         continue;
            %     end
            %
            %     --- C. 生成独立的目标 T3 (为了满足增益比 < 5) ---
            %     min_allowed_R_for_T3 = max(R1, R2) / sqrt(max_amp_ratio);
            %     max_allowed_R_for_T3 = min(R1, R2) * sqrt(max_amp_ratio);
            %
            %     R3_min_bound = max([2 * distance_res, min_allowed_R_for_T3]);
            %     R3_max_bound = min([R_max - distance_res, max_allowed_R_for_T3]);
            %
            %     if R3_min_bound >= R3_max_bound
            %         continue; % 无法同时满足边界条件，重新生成
            %     end
            %
            %     R3 = R3_min_bound + rand * (R3_max_bound - R3_min_bound);
            %
            %     确保 T3 在距离上与 T1, T2 至少间隔 1 个分辨单元
            %     if abs(R3 - R1) < 1.0 * distance_res || abs(R3 - R2) < 1.0 * distance_res
            %         continue;
            %     end
            %
            %     V3 = (rand * 2 - 1) * V_max * 0.9;
            %     确保 T3 在速度上与 T1, T2 至少间隔 1 个分辨单元
            %     if abs(V3 - V1) < 1.0 * v_res || abs(V3 - V2) < 1.0 * v_res
            %         continue;
            %     end
            %
            %     整合并排序
            %     T_Range = [R1, R2, R3];
            %     T_Velocity = [V1, V2, V3];
            %
            %     [T_Range, sort_idx] = sort(T_Range);
            %     T_Velocity = T_Velocity(sort_idx);

            % --- D. 计算幅值并做最终验证 ---
            raw_trend = 1 ./ (T_Range.^2);
            current_power = mean(raw_trend.^2);
            true_h_abs = raw_trend / sqrt(current_power); % 归一化

            min_h = min(true_h_abs);
            max_h = max(true_h_abs);
            amp_ratio = max_h / min_h;

            % 【核心条件验证】 增益比 < 5
            if min_h > 0.0065 && amp_ratio < max_amp_ratio
                is_valid = true;
                true_h = true_h_abs .* exp(1j * 2 * pi * rand(1, num_targets));
            end

            if loop_safety_count > 50000
                warning('无法生成满足条件的目标，已跳出循环。');
                break;
            end
        end

        true_li = T_Range ./ distance_res;
        true_ki = T_Velocity ./ v_res;

        %% === 2. 生成发射信号 dd (动态功率归一化版) ===
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
                dd(row, (j-1)*g + active_idx) = sym_val; % 去掉硬编码的 factor0
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
                    dd(row, guard_cols_range(i_col + act_idx - 1)) = sym_val; % 去掉硬编码的 factor0
                end
                i_col = i_col + grp_sz;
            end
        end

        % --- 【关键修改】动态功率归一化 (根据当前 SNR) ---
        current_signal_power = sum(abs(dd(:)).^2) / (M * N);
        target_power = 10^(current_SNR / 10);
        dd = dd * sqrt(target_power / current_signal_power);

        %% === 3. 生成接收信号 (Ground Truth) ===
        Y_clean = zeros(M, N);
        for p = 1:num_targets
            Y_clean = Y_clean + generate_echo_physical_truth(dd, true_h(p), true_li(p), true_ki(p), M, N);
        end
        noise = sqrt(1/2) * (randn(M, N) + 1j * randn(M, N));
        Y = Y_clean + noise;

        %% === 4. MTASA 检测算法 (修正后) ===
        [est_li, est_ki, est_h] = run_mtasa_detection(Y, dd, M, N, num_targets, delta_f);

        %% === 5. 误差计算与目标配对 ===
        est_dist = est_li * distance_res;
        est_vel  = est_ki * v_res;

        % 1. 构建误差矩阵
        dist_err_mat = zeros(num_targets, num_targets);
        vel_err_mat  = zeros(num_targets, num_targets);
        for t = 1:num_targets
            for e = 1:num_targets
                % 距离误差 (考虑周期性)
                l_true = true_li(t); l_est = est_li(e);
                diff_l = l_est - l_true;
                diff_l = diff_l - M * round(diff_l / M);
                dist_err_mat(t, e) = abs(diff_l * distance_res);

                % 速度误差 (考虑周期性)
                k_true = true_ki(t); k_est = est_ki(e);
                diff_k = k_est - k_true;
                diff_k = diff_k - N * round(diff_k / N);
                vel_err_mat(t, e) = abs(diff_k * v_res);
            end
        end

        % 2. 应用门限 (Gating)
        valid_threshold_dist = 1.0 * distance_res; % 距离门限：1个分辨单元
        cost_mat = dist_err_mat;
        cost_mat(dist_err_mat > valid_threshold_dist) = inf;

        % 3. 贪心选择 (Greedy Selection)
        matched_pairs = []; % [true_idx, est_idx]
        while true
            [min_val, min_linear_idx] = min(cost_mat(:));
            if isinf(min_val)
                break; % 剩下的匹配不上或已匹配完
            end
            [row, col] = ind2sub(size(cost_mat), min_linear_idx);
            matched_pairs = [matched_pairs; row, col];
            cost_mat(row, :) = inf;
            cost_mat(:, col) = inf;
        end

        % 4. 统计本次 MC 的误差
        mc_range_sq_err = 0;
        mc_vel_sq_err = 0;
        mc_valid_count = size(matched_pairs, 1);

        for i = 1:mc_valid_count
            t_idx = matched_pairs(i, 1);
            e_idx = matched_pairs(i, 2);
            r_err = dist_err_mat(t_idx, e_idx);
            v_err = vel_err_mat(t_idx, e_idx);
            mc_range_sq_err = mc_range_sq_err + r_err^2;
            mc_vel_sq_err   = mc_vel_sq_err + v_err^2;
        end

        % 累加到总统计
        total_valid_targets_SNR = total_valid_targets_SNR + mc_valid_count;
        sum_sq_err_range = sum_sq_err_range + mc_range_sq_err;
        sum_sq_err_vel   = sum_sq_err_vel + mc_vel_sq_err;
        dist(mc) = mc_range_sq_err;

        %% === 【关键修改】修复打印输出错位 ===
        fprintf('\n--- 目标匹配与估计结果 ---\n');
        fprintf('%-11s | %-27s | %-27s | %-27s\n', '配对状态', '距离维度 (Range Index)', '速度维度 (Velocity Index)', '增益幅度 (Gain Amplitude)');
        fprintf('%-11s | %-8s %-8s %-8s | %-8s %-8s %-8s | %-8s %-8s %-8s\n', ...
            'T <-> E', '真实值', '估计值', '偏差', '真实值', '估计值', '偏差', '真实值', '估计值', '偏差');

        % 1. 找出未成功配对的真实目标和估计目标
        all_idx = 1:num_targets;
        if ~isempty(matched_pairs)
            unmatched_T = setdiff(all_idx, matched_pairs(:, 1)');
            unmatched_E = setdiff(all_idx, matched_pairs(:, 2)');
        else
            unmatched_T = all_idx;
            unmatched_E = all_idx;
        end

        % 2. 对未配对的目标进行“强制就近配对”（仅用于打印展示偏差，不计入正确统计）
        forced_pairs = [];
        if ~isempty(unmatched_T) && ~isempty(unmatched_E)
            temp_cost = dist_err_mat(unmatched_T, unmatched_E); % 调取原始的真实误差矩阵
            while true
                [min_val, min_linear_idx] = min(temp_cost(:));
                if isinf(min_val) || isempty(min_val)
                    break;
                end
                [r, c] = ind2sub(size(temp_cost), min_linear_idx);
                forced_pairs = [forced_pairs; unmatched_T(r), unmatched_E(c)];
                temp_cost(r, :) = inf;
                temp_cost(:, c) = inf;
            end
        end

        % 3. 合并成功配对和强制配对，并打上标记状态 (1=成功, 0=失败/超限)
        all_display_pairs = [];
        if ~isempty(matched_pairs)
            all_display_pairs = [matched_pairs, ones(size(matched_pairs, 1), 1)];
        end
        if ~isempty(forced_pairs)
            all_display_pairs = [all_display_pairs; forced_pairs, zeros(size(forced_pairs, 1), 1)];
        end

        % 按真实目标(T)序号升序排序，方便查看
        if ~isempty(all_display_pairs)
            all_display_pairs = sortrows(all_display_pairs, 1);
        end

        % 4. 统一打印所有结果
        for i = 1:size(all_display_pairs, 1)
            t_idx = all_display_pairs(i, 1);
            e_idx = all_display_pairs(i, 2);
            is_success = all_display_pairs(i, 3);

            % 获取数据
            r_true = true_li(t_idx); r_est = est_li(e_idx); r_diff = r_est - r_true;
            v_true = true_ki(t_idx); v_est = est_ki(e_idx); v_diff = v_est - v_true;
            h_true_amp = abs(true_h(t_idx)); h_est_amp = abs(est_h(e_idx)); h_diff = h_est_amp - h_true_amp;

            % 正常配对用 <->，超出误差门限的出问题配对用 <!> 标记
            if is_success
                pair_str = sprintf('T[%d]<->E[%d]', t_idx, e_idx);
            else
                pair_str = sprintf('T[%d]<!>E[%d]', t_idx, e_idx);
            end

            fprintf('%-11s | %8.5f %8.5f %+8.5f | %8.5f %8.5f %+8.5f | %8.5f %8.5f %+8.5f\n', ...
                pair_str, r_true, r_est, r_diff, v_true, v_est, v_diff, h_true_amp, h_est_amp, h_diff);
        end

        if size(matched_pairs, 1) < num_targets
            warning('本次实验有目标未能成功配对（超出误差门限）！带 "<!>" 标记的为就近强制配对结果，供排查极大误差参考。');
        end

        % 绘图 (热力图) 建议调试完毕后将这部分注释掉以加快仿真速度
        figure('Name', 'DD Domain Range-Doppler Heatmap', 'Color', 'w');
        imagesc(1:N, 1:M, abs(Y));
        colormap('jet'); colorbar;
        xlabel('Doppler Index (Column Index)');
        ylabel('Delay Index (Row Index)');
        title(['Received Signal Y (Magnitude) | SNR=' num2str(current_SNR) 'dB']);
        hold on;
        for t = 1:num_targets
            y_pos = true_li(t) + 1;
            k_val = true_ki(t);
            if k_val < 0
                x_pos = N + k_val + 1;
            else
                x_pos = k_val + 1;
            end
            plot(x_pos, y_pos, 'rx', 'MarkerSize', 12, 'LineWidth', 2.5);
            text(x_pos, y_pos - 1, sprintf(' T%d', t), 'Color', 'white', 'FontWeight', 'bold');
        end
        hold off;
        set(gca, 'YDir', 'reverse');
        fprintf('MC=%d 时的真实相位: T1=%.4f, T2=%.4f, T3=%.4f (弧度)\n', mc, angle(true_h(1)), angle(true_h(2)), angle(true_h(3)));
    end

    %% === 恢复计算检测概率与 MSE ===
    if total_valid_targets_SNR > 0
        MSE_Range(s_idx) = sqrt(sum_sq_err_range / total_valid_targets_SNR);
        MSE_Velocity(s_idx) = sqrt(sum_sq_err_vel / total_valid_targets_SNR);
    else
        MSE_Range(s_idx) = NaN;
        MSE_Velocity(s_idx) = NaN;
    end

    Detection_Prob(s_idx) = total_valid_targets_SNR / (num_monte_carlo * num_targets);

    % fprintf(' => 当前SNR完成! 有效匹配数: %d, RMSE_Range: %.4e, 检测概率: %.2f%%\n', ...
    %     total_valid_targets_SNR, MSE_Range(s_idx), Detection_Prob(s_idx)*100);


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
function [est_li, est_ki, est_h] = run_mtasa_detection(Y, dd, M, N, num_targets, delta_f)
% MTASA: Multi-Target Active Sensing Algorithm
est_li = zeros(1, num_targets);
est_ki = zeros(1, num_targets);
est_h  = zeros(1, num_targets);

% --- [SIC] 初始化 ---
Residual = Y;
for k = 1:num_targets
    [h_hat, l_hat, k_hat] = two_stage_search_matched_golden(Residual, dd, M, N);
    est_li(k) = l_hat;
    est_ki(k) = k_hat;
    est_h(k)  = h_hat;
    Signal_Est = generate_echo_physical_truth(dd, h_hat, l_hat, k_hat, M, N);
    Residual = Residual - Signal_Est;
end

% --- [PIC] 并行干扰消除迭代 ---
MaxIter = 100;
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

% 后处理
mask_neg = est_ki > (N/2);
est_ki(mask_neg) = est_ki(mask_neg) - N;
est_li = mod(est_li, M);
end

function [h_hat, l_hat, k_hat] = two_stage_search_matched_golden(Y_received, dd, M, N)
% 粗搜索
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

% 精搜索 - 二维黄金分割法
delta_f = 120e3;
T = 1 / delta_f;
tau_c = l_coarse_idx / (M * delta_f);
nu_c  = k_coarse_idx / (N * T);

res_tau = 1 / (M * delta_f);
res_nu  = 1 / (N * T);

% --- 【关键修改】区间收紧为 0.5 ---
a_l = tau_c - 1 * res_tau; a_u = tau_c +1 * res_tau;
b_l = nu_c - 1 * res_nu;   b_u = nu_c + 1 * res_nu;

mu = (sqrt(5) - 1) / 2;
Iter = 20;
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