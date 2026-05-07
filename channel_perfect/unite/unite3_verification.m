%% 剥离验证实验：已知真实索引下的纯符号解调性能对比 (采用安全 waitbar 进度条)
clc; clear; close all;
M=64;
N=32;
fc=64e9; delta_f=120e3;
u_max=120*1000/3600; c=3e8;
% 计算抽头
delay=[30,150,310,370,710,1090,1730,2510]*1e-9;
delay_tap0=delay*M*delta_f;
li=round(delay_tap0);
doppler_max=fc*u_max/c;
k_max0=doppler_max*N*(1/delta_f);

% EVA信道，时延为0的路径是增益最大的，以后依次递减
h_p_db=[1.5,1.4,3.6,0.6,9.1,7,12,16.9];
relative_power_linear = 10.^((-1) * h_p_db/10); % db转换为功率
total_power = sum(relative_power_linear);
h_p = relative_power_linear / total_power; % 每一个增益所占据的百分比
l_max=20; k_max=2;
P=length(delay); % 路径数
cen_modu = 2;   % BPSK
guard_modu = 4; % QPSK
bits_per_qam_symbol = log2(cen_modu);
bits_per_guard_symbol=log2(guard_modu);
g = 4; % IM分组参数
guard_data_power_factor = 1;

% 计算容量与激活符号数
mo=0;
q1 = floor((N-4*k_max-1)/ g);
bit_row=(log2(guard_modu)+mo)+q1*(log2(guard_modu)+log2(nchoosek(g,1))); % guard一行bit数
bit_row_num=(1+q1)*(2*l_max+1); % guard非零dd块个数
cen_num=(N/g)*(M-(1+2*l_max));
ans3=cen_num+bit_row_num;
bit=(log2(cen_modu)+log2(nchoosek(g,1)))*(cen_num)+(2*l_max+1)*bit_row;

% 仿真参数设定
SNR=10:2:18;
% 为了快速验证趋势，迭代次数稍微调小，若要跑平滑曲线可调回您的原参数
iter=[4e4, 4e4, 4e5, 4e5, 4e6];
power_persym=10.^(SNR/10); % 每个符号的功率
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
power=power_persym*ans2; % 基准总功率
factor=power/ans3;
xp=sqrt(power_persym(1)*1e5); % 导频能量

% 用于保存两种极端条件下的误码率
ber_perf = zeros(1, length(SNR)); % 完美信道 + 纯符号解调
ber_imp = zeros(1, length(SNR));  % 非完美信道 + 纯符号解调

% ================= 新增：安全进度条初始化 =================
total_frames = sum(ceil(iter ./ bit));
current_frame = 0;
h_wait = waitbar(0, '准备开始仿真...', 'Name', '剥离验证实验进度');
% =========================================================

for i_snr=1:length(SNR)
    num_frames = ceil(iter(i_snr)/bit);
    
    frame_err_perf = zeros(1, num_frames);
    frame_err_imp  = zeros(1, num_frames);
    frame_bits = zeros(1, num_frames);
    
    %% 功率分配
    data_zone_mask = false(M,N);
    data_rows = (l_max+2):(M-l_max); % 数据行号
    data_zone_mask(data_rows, :) = true;
    
    num_central_symbols = cen_num;
    num_guard_symbols = ans3-num_central_symbols;
    
    denominator = num_central_symbols + num_guard_symbols * guard_data_power_factor;
    power_per_central_symbol = power(i_snr) / denominator;
    power_per_guard_symbol = power_per_central_symbol * guard_data_power_factor;
    
    scale_central = sqrt(power_per_central_symbol);
    scale_guard = sqrt(power_per_guard_symbol);
    
    for f=1:num_frames
        
        % 动态刷新进度条
        current_frame = current_frame + 1;
        if mod(f, max(1, floor(num_frames/20))) == 0 || f == num_frames
            wait_msg = sprintf('正在运行 SNR = %d dB: 第 %d / %d 帧', SNR(i_snr), f, num_frames);
            waitbar(current_frame / total_frames, h_wait, wait_msg);
        end
        
        % 随机信道参数
        theta0= -pi + 2*pi*rand(1, length(delay));
        doppler=doppler_max*cos(theta0) ;
        ki=doppler*N*(1/delta_f);
        h_exp=exp( 1j*2 * pi * rand(1, length(h_p_db)));
        
        %% 发射端：基于信息比特生成OTFS-IM信号
        dd = zeros(M, N);
        original_bits_stream = [];
        num_groups_data = N / g; 
        bits_per_IM_data = log2(nchoosek(g, 1));
        
        % --- A: 处理中央数据区域 ---
        for row = data_rows
            for j = 1:num_groups_data
                im_bits = randi([0 1], 1, bits_per_IM_data);
                qam_bits = randi([0 1], 1, bits_per_qam_symbol);
                original_bits_stream = [original_bits_stream, im_bits, qam_bits];
                
                active_local_idx = bi2de(im_bits, 'left-msb') + 1;
                symbol_int = bi2de(qam_bits, 'left-msb');
                qam_symbol = qammod(symbol_int, cen_modu, 'UnitAveragePower', true);
                
                start_col = (j - 1) * g + 1; 
                dd(row, start_col + active_local_idx - 1) = qam_symbol;
            end
        end
        
        % --- B: 处理保护间隔区域 ---
        guard_rows = [1:(l_max + 1), (M - l_max + 1):M];
        guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
        num_guard_cols = length(guard_cols_range);
        for row = guard_rows
            i_col = 1; 
            while i_col <= num_guard_cols 
                current_group_size = min(g, num_guard_cols - i_col + 1);
                if current_group_size > 1 
                    if current_group_size == g
                        num_im_bits = floor(log2(nchoosek(g, 1)));
                        im_bits = randi([0 1], 1, num_im_bits);
                        active_local_idx = bi2de(im_bits, 'left-msb') + 1;
                    else
                        active_local_idx = randi(current_group_size);
                        im_bits = zeros(1, current_group_size);
                        im_bits(active_local_idx) = 1;
                    end
                    qam_bits = randi([0 1], 1, bits_per_guard_symbol);
                    original_bits_stream = [original_bits_stream, im_bits, qam_bits];
                    
                    symbol_int = bi2de(qam_bits, 'left-msb');
                    qam_symbol = qammod(symbol_int, guard_modu, 'UnitAveragePower', true);
                    global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
                    dd(row, global_col_indices(active_local_idx)) = qam_symbol;
                end
                i_col = i_col + current_group_size;
            end
        end
        
        % 幅度缩放与导频插入
        dd(data_rows,:) = dd(data_rows,:)* scale_central;
        guard_data_indices=[1:(l_max+1), (M-l_max+1):M];
        dd(guard_data_indices,:) = dd(guard_data_indices,:)* scale_guard;
        dd(1,1)=xp;
        
        %% 信道传输
        hw=zeros(M,N); 
        for l=0:M-1
            for k=0:N-1
                for i=1:P
                    theta=exp(1j*2*pi*ki(i)*(li(i))/(M*N));
                    delta_term = (l == li(i)); 
                    hw(l+1,k+1)=hw(l+1,k+1)+h_p(i)*h_exp(i)*delta_term*zeta_N(k-ki(i),N).*theta;
                end
            end
        end
        H_freq = fft2(hw);
        dd_freq = fft2(dd);
        y_freq = H_freq .* dd_freq;
        y = ifft2(y_freq);
        
        % 加噪
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));
        Y=y+reshape(noise,M,N);
        
        %% (1) 信道估计 (粗搜+精搜)
        YABS=abs(Y);
        [max_values0, ~] = max(YABS, [], 2);
        
        step=0.1;
        colnum=length(-k_max:step:k_max);
        r=zeros(l_max+1,colnum);
        for l=1:l_max+1
            jj=1; 
            yd=Y(l,:);
            for stepval=-k_max:step:k_max  
                for k=1:size(Y,2)
                    r(l,jj)=r(l,jj)+yd(k)*conj(zeta_N(k-stepval-1,N));
                end
                jj=jj+1; 
            end
        end
        
        [max_values, max_indices] = max(abs(r), [], 2);
        [~, sorted_indices] = sort(max_values0(1:l_max+1), 'descend');
        top_indices = sorted_indices(1:8);
        li_est = sort(top_indices - 1);
        h_est=[]; h_phi_est=[]; ki_est=[];
        for i=1:length(li_est)
            yd= Y(li_est(i)+1,:); 
            peak_idx = max_indices(li_est(i)+1);
            left_idx = max(1, peak_idx - 1);
            right_idx = min(colnum, peak_idx + 1);
            b_l = (-k_max) + step * (left_idx - 1);
            b_u = (-k_max) + step * (right_idx - 1);
            
            % GSS 极简版
            max_iter=30;
            eta = (sqrt(5)-1)/2;
            b1 = b_u - eta * (b_u - b_l); 
            b2 = b_l + eta * (b_u - b_l);
            
            f1_val = 0; f2_val = 0;
            for k_idx = 1:N
                f1_val = f1_val + yd(k_idx) * conj(zeta_N(k_idx-b1-1, N));
                f2_val = f2_val + yd(k_idx) * conj(zeta_N(k_idx-b2-1, N));
            end
            f1_val = abs(f1_val); 
            f2_val = abs(f2_val);
            for iter_gss = 1:max_iter
                if f1_val > f2_val
                    b_u = b2; b2 = b1; f2_val = f1_val;
                    b1 = b_u - eta * (b_u - b_l);
                    f1_val = 0;
                    for k_idx = 1:N
                        f1_val = f1_val + yd(k_idx) * conj(zeta_N(k_idx-b1-1, N));
                    end
                    f1_val = abs(f1_val);
                else
                    b_l = b1; b1 = b2; f1_val = f2_val;
                    b2 = b_l + eta * (b_u - b_l);
                    f2_val = 0;
                    for k_idx = 1:N
                        f2_val = f2_val + yd(k_idx) * conj(zeta_N(k_idx-b2-1, N));
                    end
                    f2_val = abs(f2_val);
                end
            end
            nu_opt = (b_l + b_u) / 2;
            ki_est=[ki_est,nu_opt];
            
            r_refined = 0;
            for k_idx = 1:N
                r_refined = r_refined + yd(k_idx) * conj(zeta_N(k_idx-nu_opt-1, N));
            end
            
            h_est = [h_est, abs(r_refined) / xp];
            if abs(r_refined) > 1e-9
                h_phi_est = [h_phi_est, r_refined / abs(r_refined)]; 
            else
                h_phi_est = [h_phi_est, 1]; 
            end
        end
        
        hw_est=zeros(M,N);
        for l=0:M-1
            for k=0:N-1
                for i=1:length(li_est)
                    theta=exp(1j*2*pi*ki_est(i)*(li_est(i))/(M*N));
                    delta_term = (l == li_est(i)); 
                    hw_est(l+1,k+1)=hw_est(l+1,k+1)+h_est(i)*h_phi_est(i)*delta_term*zeta_N(k-ki_est(i),N).*theta;
                end
            end
        end
        %% (2) 双路 MMSE 均衡 (完美信道 vs 非完美信道)
        alpha = 1 / factor(i_snr);
        y_freq_rx = fft2(Y); 
        
        % 路径 1: 非完美信道估计的均衡
        H_freq_imp = fft2(hw_est);
        H_mmse_freq_imp = conj(H_freq_imp) ./ (abs(H_freq_imp).^2 + alpha);
        dd_est_matrix_imp = ifft2(H_mmse_freq_imp .* y_freq_rx);
        
        % 路径 2: 完美信道的均衡
        H_freq_perf = fft2(hw);
        H_mmse_freq_perf = conj(H_freq_perf) ./ (abs(H_freq_perf).^2 + alpha);
        dd_est_matrix_perf = ifft2(H_mmse_freq_perf .* y_freq_rx);
        
        %% (3) 剥离验证：上帝视角解调 (已知真实激活位置)
        est_bits_imp = [];
        est_bits_perf = [];
        
        % --- A: 中央区域 ---
        for row = data_rows
            for j = 1:num_groups_data
                start_col = (j - 1) * g + 1;
                group_cols = start_col : start_col + g - 1;
                
                % 获取真实的发送组，找到真实激活索引
                true_group_tx = dd(row, group_cols);
                true_local_idx = find(true_group_tx ~= 0);
                if isempty(true_local_idx); true_local_idx = 1; end % 防错
                true_local_idx = true_local_idx(1); 
                
                % 提取该位置的信号
                val_imp  = dd_est_matrix_imp(row, group_cols(true_local_idx)) / scale_central;
                val_perf = dd_est_matrix_perf(row, group_cols(true_local_idx)) / scale_central;
                
                % 纯符号硬判决
                det_sym_int_imp  = qamdemod(val_imp, cen_modu, 'UnitAveragePower', true, 'OutputType', 'integer');
                det_sym_int_perf = qamdemod(val_perf, cen_modu, 'UnitAveragePower', true, 'OutputType', 'integer');
                
                % 组合比特
                im_bits_true = de2bi(true_local_idx - 1, bits_per_IM_data, 'left-msb');
                qam_bits_imp  = de2bi(det_sym_int_imp, bits_per_qam_symbol, 'left-msb');
                qam_bits_perf = de2bi(det_sym_int_perf, bits_per_qam_symbol, 'left-msb');
                
                est_bits_imp  = [est_bits_imp, im_bits_true, qam_bits_imp];
                est_bits_perf = [est_bits_perf, im_bits_true, qam_bits_perf];
            end
        end
        
        % --- B: 保护间隔区域 ---
        for row = guard_rows
            i_col = 1;
            while i_col <= num_guard_cols
                current_group_size = min(g, num_guard_cols - i_col + 1);
                if current_group_size  > 1
                    global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
                    
                    % 获取真实的激活位置
                    true_group_tx = dd(row, global_col_indices);
                    true_local_idx = find(true_group_tx ~= 0);
                    if isempty(true_local_idx); true_local_idx = 1; end
                    true_local_idx = true_local_idx(1);
                    
                    % 提取该位置信号
                    val_imp  = dd_est_matrix_imp(row, global_col_indices(true_local_idx)) / scale_guard;
                    val_perf = dd_est_matrix_perf(row, global_col_indices(true_local_idx)) / scale_guard;
                    
                    % 纯符号硬判决
                    det_sym_int_imp  = qamdemod(val_imp, guard_modu, 'UnitAveragePower', true, 'OutputType', 'integer');
                    det_sym_int_perf = qamdemod(val_perf, guard_modu, 'UnitAveragePower', true, 'OutputType', 'integer');
                    
                    % 组合比特
                    if current_group_size == g
                        im_bits_true = de2bi(true_local_idx - 1, bits_per_IM_data, 'left-msb');
                    else
                        im_bits_true = zeros(1, current_group_size);
                        im_bits_true(true_local_idx) = 1;
                    end
                    qam_bits_imp  = de2bi(det_sym_int_imp, bits_per_guard_symbol, 'left-msb');
                    qam_bits_perf = de2bi(det_sym_int_perf, bits_per_guard_symbol, 'left-msb');
                    
                    est_bits_imp  = [est_bits_imp, im_bits_true, qam_bits_imp];
                    est_bits_perf = [est_bits_perf, im_bits_true, qam_bits_perf];
                end
                i_col = i_col + current_group_size;
            end
        end
        
        % 计算这一帧的错误数
        [err_imp, ~] = biterr(original_bits_stream, est_bits_imp);
        [err_perf, ~] = biterr(original_bits_stream, est_bits_perf);
        
        frame_err_imp(f)  = err_imp;
        frame_err_perf(f) = err_perf;
        frame_bits(f) = length(original_bits_stream);
    end
    
    % 统计该SNR下的平均误码率
    ber_imp(i_snr)  = sum(frame_err_imp) / sum(frame_bits);
    ber_perf(i_snr) = sum(frame_err_perf) / sum(frame_bits);
    
    fprintf('>>> SNR = %d 完成 | BER_Imp: %e | BER_Perf: %e\n', SNR(i_snr), ber_imp(i_snr), ber_perf(i_snr));
end

% 关闭进度条
if isvalid(h_wait)
    close(h_wait);
end

%% 绘制剥离实验结果图
figure('Color', 'w');
semilogy(SNR, ber_perf, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(SNR, ber_imp, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('信噪比 (dB)', 'FontSize', 12);
ylabel('纯符号比特误码率 (BER)', 'FontSize', 12);
title('剥离验证：已知真实索引下的纯符号解调性能', 'FontSize', 14);
legend('完美信道估计 (已知索引)', '非完美信道估计 (已知索引)', 'Location', 'SouthWest');

% -------------------------------------------------------------
% 【修复版】zeta_N 函数 - 支持全向量化，避免标量 || 报错
% -------------------------------------------------------------
function output = zeta_N(k, N)
    tolerance = 1e-10; 
    output = zeros(size(k));
    
    % 使用逻辑索引，支持数组级判断
    idx_zero = (abs(k) < tolerance) | (abs(k - N) < tolerance);
    output(idx_zero) = 1;
    
    idx_norm = ~idx_zero;
    if any(idx_norm(:))
        k_norm = k(idx_norm);
        exp_term = exp(-1i * pi * (N - 1) * k_norm / N);
        numerator = sin(pi * k_norm);
        numerator(abs(numerator) < tolerance) = 0; % 过滤极小浮点误差
        denominator = sin(pi * k_norm / N);
        output(idx_norm) = exp_term .* (numerator ./ (N * denominator));
    end
end