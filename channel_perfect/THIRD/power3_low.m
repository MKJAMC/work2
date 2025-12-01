%% 索引比特的进阶检测方法 (O(g) 复杂度优化版)
clc;clear;close all;
M=64;
N=32;
fc=64e9;delta_f=120e3;
u_max=120*1000/3600;c=3e8;
%计算抽头
delay=[30,150,310,370,710,1090,1730,2510]*1e-9;
delay_tap0=delay*M*delta_f;
li=round(delay_tap0);
doppler_max=fc*u_max/c;
k_max0=doppler_max*N*(1/delta_f);
% EVA信道
h_p_db=[1.5,1.4,3.6,0.6,9.1,7,12,16.9];
relative_power_linear = 10.^((-1) * h_p_db/10);%db转换为功率
total_power = sum(relative_power_linear);
h_p = relative_power_linear / total_power;%每一个增益所占据的百分比
l_max=20;k_max=2;
P=length(delay);%路径数
cen_modu =2;
guard_modu= 4;
bits_per_qam_symbol = log2(cen_modu);
bits_per_guard_symbol=log2(guard_modu);
g = 4; % IM分组参数
guard_data_power_factor =1;
% 速率=位置信息+调制阶数
mo=mod((N-4*k_max-1),g);
q1 = floor((N-4*k_max-1)/ g);
%一行的速率与绿色区间激活个数
bit_row=(log2(guard_modu)+mo)+q1*(log2(guard_modu)+log2(nchoosek(g,1)));%guard一行bit数
bit_row_num=(1+q1)*(2*l_max+1);%guard非零dd块个数
cen_num=(N/g)*(M-(1+2*l_max));
ans3=cen_num+bit_row_num;
bit=(log2(cen_modu)+log2(nchoosek(g,1)))*(cen_num)+(2*l_max+1)*bit_row;
% 符号检测SNR
SNR=10:2:16;
iter=[4e4,4e4,4e5,4e5];
% SNR=18;iter=[4e6];
cen_constellation = qammod(0:cen_modu-1, cen_modu, 'UnitAveragePower', true);
guard_constellation = qammod(0:guard_modu-1, guard_modu, 'UnitAveragePower', true);

%% 主循环
power_persym=10.^(SNR/10);%每个符号的功率
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
power=power_persym*ans2;%ans2为基准，求得总功率
ber = zeros(1, length(SNR));
factor=power/ans3;

for i_snr=1:length(SNR)
 
    ber_iter = zeros(1, iter(i_snr));
    num_frames =ceil(iter(i_snr)/bit);
    frame_errors = zeros(1, num_frames);
    frame_bits = zeros(1, num_frames);
    
    %% 功率分配
    data_rows = (l_max+2):(M-l_max); 
    num_central_symbols = cen_num;
    num_guard_symbols = ans3-num_central_symbols;
    denominator = num_central_symbols + num_guard_symbols * guard_data_power_factor;
    if denominator > 0
        power_per_central_symbol = power(i_snr) / denominator;
    else
        power_per_central_symbol = 0; 
    end
    power_per_guard_symbol = power_per_central_symbol * guard_data_power_factor;
    scale_central = sqrt(power_per_central_symbol);
    scale_guard = sqrt(power_per_guard_symbol);

    for f=1:num_frames
        disp(f) % 减少打印以提高速度
        theta0= -pi + 2*pi*rand(1, length(delay));
        doppler=doppler_max*cos(theta0) ;
        ki=doppler*N*(1/delta_f);
        h_exp=exp( 1j*2 * pi * rand(1, length(h_p_db)));
        
        %% 发射端
        dd = zeros(M, N);
        original_bits_stream = [];
        data_rows = (l_max + 2):(M - l_max);
        num_groups_data = N / g; 
        bits_per_IM_data = log2(nchoosek(g, 1));
        
        % --- A: 发射端 中央数据区域 ---
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
        
        % --- B: 发射端 保护间隔区域 ---
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

        dd(data_rows,:) = dd(data_rows,:)* scale_central;
        guard_data_indices=[1:(l_max+1), (M-l_max+1):M];
        dd(guard_data_indices,:) = dd(guard_data_indices,:)* scale_guard;
        
        % 信道
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
        
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));
        y=y+reshape(noise,M,N);
        
        %% 低复杂度MMSE与全搜索检测
        y_freq = fft2(y); 
        alpha = 1 / factor(i_snr);
        H_mmse_freq = conj(H_freq) ./ (abs(H_freq).^2 + alpha);
        z_freq_initial = H_mmse_freq .* y_freq;
        dd_est_matrix = ifft2(z_freq_initial);
        
        estimated_non_zero_indices=[];estimated_bits_stream=[];
        
        % ============================================================
        % 优化区域 A: 中央数据区 - O(g) 复杂度全搜索
        % ============================================================
        for row = data_rows
            for j = 1:num_groups_data
                start_col = (j - 1) * g + 1;
                group_cols = start_col : start_col + g - 1;
                group_est = dd_est_matrix(row, group_cols) / scale_central;
                
                % [优化1] 预先计算当前分组的总能量 (O(g))
                group_energy = sum(abs(group_est).^2);
                
                max_metric = -inf;
                best_local_idx = -1;
                best_symbol = 0;
                
                % 遍历所有可能激活的位置 k (O(g))
                for k = 1:g
                    % [优化2] 直接硬判决得到当前位置的最优符号 (O(1))
                    % 替代了原本的星座图遍历 O(M)
                    detected_symbol_int = qamdemod(group_est(k), cen_modu, 'UnitAveragePower', true, 'OutputType', 'integer');
                    current_symbol = qammod(detected_symbol_int, cen_modu, 'UnitAveragePower', true);
                    
                    % [优化3] 利用预计算的总能量计算度量 (O(1))
                    % Metric = -|y_k - s|^2 - sum(|y_m|^2)
                    %        = -|y_k - s|^2 - (Total_Energy - |y_k|^2)
                    match_metric = -abs(group_est(k) - current_symbol)^2;
                    penalty_metric = group_energy - abs(group_est(k))^2;
                    
                    metric = match_metric - penalty_metric;
                    
                    if metric > max_metric
                        max_metric = metric;
                        best_local_idx = k;
                        best_symbol = current_symbol;
                    end
                end
                
                % 解码结果
                detected_local_idx = best_local_idx;
                detected_symbol = best_symbol;
                
                im_bits_est = de2bi(detected_local_idx - 1, bits_per_IM_data, 'left-msb');
                qam_bits_est = qamdemod(detected_symbol, cen_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
                estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
            end
        end
        
        % ============================================================
        % 优化区域 B: 保护间隔区 - O(g) 复杂度全搜索
        % ============================================================
        for row = guard_rows
            i_col = 1;
            while i_col <= num_guard_cols
                current_group_size = min(g, num_guard_cols - i_col + 1);
                if current_group_size  > 1
                    global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
                    group_est = dd_est_matrix(row, global_col_indices) / scale_guard;
                    
                    % [优化1] 预计算总能量
                    group_energy = sum(abs(group_est).^2);
                    
                    max_metric = -inf;
                    best_local_idx = -1;
                    best_symbol = 0;
                    
                    % 遍历位置
                    for k = 1:current_group_size
                        % [优化2] 硬判决
                        detected_symbol_int = qamdemod(group_est(k), guard_modu, 'UnitAveragePower', true, 'OutputType', 'integer');
                        current_symbol = qammod(detected_symbol_int, guard_modu, 'UnitAveragePower', true);
                        
                        % [优化3] 计算度量
                        match_metric = -abs(group_est(k) - current_symbol)^2;
                        penalty_metric = group_energy - abs(group_est(k))^2;
                        
                        metric = match_metric - penalty_metric;
                        
                        if metric > max_metric
                            max_metric = metric;
                            best_local_idx = k;
                            best_symbol = current_symbol;
                        end
                    end
                    
                    detected_local_idx = best_local_idx;
                    detected_symbol = best_symbol;
                    
                    % 解码
                    if current_group_size == g
                        im_bits_est = de2bi(detected_local_idx - 1, bits_per_IM_data, 'left-msb');
                    else 
                        im_bits_est = zeros(1, current_group_size);
                        im_bits_est(detected_local_idx) = 1;
                    end
                    qam_bits_est = qamdemod(detected_symbol, guard_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
                    estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
                end
                i_col = i_col + current_group_size;
            end
        end
        
        [num_errors_this_frame, ~] = biterr(original_bits_stream, estimated_bits_stream);
        num_bits_this_frame = length(original_bits_stream);
        frame_errors(f) = num_errors_this_frame;
        frame_bits(f) = num_bits_this_frame;
    end
    total_bit_errors = sum(frame_errors);
    total_bits_transmitted = sum(frame_bits);
    ber(i_snr) =total_bit_errors / total_bits_transmitted;
end

figure();
semilogy(SNR, ber, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('信噪比 (dB)');ylabel('比特误码率 (BER)');title('O(g) 复杂度联合度量检测性能');legend('IM-Proposed-O(g)');

function output = zeta_N(k, N)
tolerance = 1e-10; 
if (k == N||k==0)
    output = 1;
else
    exp_term = exp(-1i * pi * (N - 1) * k / N);
    numerator = sin(pi * k); 
    if abs(numerator) < tolerance
        numerator = 0;
    end
    denominator = sin(pi * k / N); 
    output = exp_term * (numerator / (N * denominator));
end
end