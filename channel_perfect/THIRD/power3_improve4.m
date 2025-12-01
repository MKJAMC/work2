%% 索引比特的进阶检测方法


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
% EVA信道，时延为0的路径是增益最大的，以后依次递减
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
%% 调试
% SNR=16;iter=4e6;

power_persym=10.^(SNR/10);%每个符号的功率
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
power=power_persym*ans2;%ans2为基准，求得总功率，方案二的功率为1
ber = zeros(1, length(SNR));
factor=power/ans3;
for i_snr=1:length(SNR)
    ber_iter = zeros(1, iter(i_snr));%一个frame错误bit数
    num_frames =ceil(iter(i_snr)/bit);
    frame_errors = zeros(1, num_frames);
    frame_bits = zeros(1, num_frames);
    %% 功率分配
    data_zone_mask = false(M,N);
    data_rows = (l_max+2):(M-l_max); % 数据行号
    data_zone_mask(data_rows, :) = true;
    % 找到在这个区域中被IM激活的符号:中央data
    num_central_symbols = cen_num;
    num_guard_symbols = ans3-num_central_symbols;
    % 解方程
    % 总数据功率 = 中央符号数*中央单位功率(默认为1) + 卫兵符号数*卫兵单位功率（1缩小guard_data_power_factor倍）
    denominator = num_central_symbols + num_guard_symbols * guard_data_power_factor;
    if denominator > 0
        power_per_central_symbol = power(i_snr) / denominator;
    else
        power_per_central_symbol = 0; % 防止除以0
    end
    power_per_guard_symbol = power_per_central_symbol * guard_data_power_factor;
    % 3.4 计算幅度缩放因子 (幅度是功率的平方根)
    scale_central = sqrt(power_per_central_symbol);
    scale_guard = sqrt(power_per_guard_symbol);


    for f=1:num_frames
        disp(f)
        theta0= -pi + 2*pi*rand(1, length(delay));
        doppler=doppler_max*cos(theta0) ;
        ki=doppler*N*(1/delta_f);
        h_exp=exp( 1j*2 * pi * rand(1, length(h_p_db)));
        %% 发射端：基于信息比特生成OTFS-IM信号
        dd = zeros(M, N);
        original_bits_stream = [];
        data_rows = (l_max + 2):(M - l_max);
        num_groups_data = N / g; % 计算分组数量
        bits_per_IM_data = log2(nchoosek(g, 1));%每组激活一个
        % --- A: 处理中央数据区域 (g = 组大小) ---
        for row = data_rows
            for j = 1:num_groups_data % 循环 N/g 次
                im_bits = randi([0 1], 1, bits_per_IM_data);%生成每一组的IM位置bit
                qam_bits = randi([0 1], 1, bits_per_qam_symbol);%生成一个qam符号的bit
                original_bits_stream = [original_bits_stream, im_bits, qam_bits];
                active_local_idx = bi2de(im_bits, 'left-msb') + 1;%将IM的bit转为位置
                symbol_int = bi2de(qam_bits, 'left-msb');
                qam_symbol = qammod(symbol_int, cen_modu, 'UnitAveragePower', true);
                start_col = (j - 1) * g + 1; % 每组的起始列
                dd(row, start_col + active_local_idx - 1) = qam_symbol;
            end
        end
        % --- B: 处理保护间隔区域 (使用独热编码, g = 组大小) ---
        guard_rows = [1:(l_max + 1), (M - l_max + 1):M];
        guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
        num_guard_cols = length(guard_cols_range);%需要索引调制的保护间隔部分
        % 遍历每一个指定的保护行
        for row = guard_rows
            i_col = 1; % 初始化保护列的起始计数器
            while i_col <= num_guard_cols % 计算当前组的大小，最大不超过 g
                current_group_size = min(g, num_guard_cols - i_col + 1);
                if current_group_size > 1 %最后余数为1，则不存在索引调制
                    if current_group_size == g
                        num_im_bits = floor(log2(nchoosek(g, 1)));
                        im_bits = randi([0 1], 1, num_im_bits);
                        active_local_idx = bi2de(im_bits, 'left-msb') + 1;
                    else
                        % --- 情况2: 剩余组 (current_group_size < g)，使用独热编码,在当前剩余的组内随机选择一个激活位置
                        active_local_idx = randi(current_group_size);
                        % 生成独热编码比特向量
                        im_bits = zeros(1, current_group_size);
                        im_bits(active_local_idx) = 1;
                    end
                    qam_bits = randi([0 1], 1, bits_per_guard_symbol);
                    original_bits_stream = [original_bits_stream, im_bits, qam_bits];
                    % 将QAM比特调制成一个复数符号
                    symbol_int = bi2de(qam_bits, 'left-msb');
                    qam_symbol = qammod(symbol_int, guard_modu, 'UnitAveragePower', true);
                    % 获取当前组在dd的实际列索引
                    global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
                    % 将QAM符号放置在被选中的激活位置上
                    dd(row, global_col_indices(active_local_idx)) = qam_symbol;
                end
                % 移动到下一个组的起始位置
                i_col = i_col + current_group_size;
            end
        end
        ans3 = nnz(dd);%总激活个数
        % 3.5 将不同的缩放因子应用到dd矩阵的不同区域
        % 注意：此时dd中的激活符号幅度还是1
        dd(data_rows,:) = dd(data_rows,:)* scale_central;
        guard_data_indices=[1:(l_max+1), (M-l_max+1):M];
        dd(guard_data_indices,:) = dd(guard_data_indices,:)* scale_guard;
        a2=abs(dd);
        % 分数多普勒域下，DD域等效信道代码,发送数据和接受数据都是dd域
        hw=zeros(M,N);Y=zeros(M,N);
        for l=0:M-1
            for k=0:N-1
                for i=1:P
                    theta=exp(1j*2*pi*ki(i)*(li(i))/(M*N));
                    delta_term = (l == li(i)); % 如果 l 等于 li(i)，则 delta 为 1
                    hw(l+1,k+1)=hw(l+1,k+1)+h_p(i)*h_exp(i)*delta_term*zeta_N(k-ki(i),N).*theta;
                end
            end
        end


        % 频域相乘
        H_freq = fft2(hw);
        % 步骤 B: 对输入 dd 做二维傅里叶变换
        dd_freq = fft2(dd);
        % 步骤 C: 在频域进行逐元素相乘 (这步等效于 H_eff * dd_vec)
        y_freq = H_freq .* dd_freq;
        % 步骤 D: 将结果通过二维逆傅里叶变换返回到原始域
        y = ifft2(y_freq);
        % 添加噪声
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪
        y=y+reshape(noise,M,N);

        %% 低复杂度MMSE
        y_freq = fft2(y); % y已经是矩阵形式，无需reshape
        alpha = 1 / factor(i_snr);
        H_mmse_freq = conj(H_freq) ./ (abs(H_freq).^2 + alpha);
        z_freq_initial = H_mmse_freq .* y_freq;
        dd_est_matrix = ifft2(z_freq_initial);
        a1=abs(dd_est_matrix);
        
        estimated_non_zero_indices=[];estimated_bits_stream=[];

        % --- A: 处理中央数据区域 (联合ML检测) ---
        for row = data_rows
            for j = 1:num_groups_data
                start_col = (j - 1) * g + 1;
                group_cols = start_col : start_col + g - 1;
                group_est = dd_est_matrix(row, group_cols) / scale_central;
                
                % 初始化度量和最佳候选
                max_metric = -inf;
                best_local_idx = -1;
                best_symbol = 0;
                
                % 遍历所有可能激活的位置 (Hypothesis k)
                for k = 1:g
                    % 遍历所有可能的QAM星座点 (Hypothesis s_q)
                    for q = 1:cen_modu
                        current_symbol = cen_constellation(q);                  
                        % 计算对数似然度量 (忽略常数项)
                        % Metric = - (distance_from_signal + distance_from_noise_zeros)
                        % y_k - s_q : 信号位置的距离
                        % y_m - 0   : 噪声位置的距离
                        metric = -abs(group_est(k) - current_symbol)^2;
                        for m = 1:g
                            if m ~= k
                                metric = metric - abs(group_est(m))^2;
                            end
                        end
                        if metric > max_metric
                            max_metric = metric;
                            best_local_idx = k;
                            best_symbol = current_symbol;
                        end
                    end
                end
                
                % 根据找到的最佳位置和符号来解码比特
                detected_local_idx = best_local_idx;
                detected_symbol = best_symbol;

                global_col_idx = start_col + detected_local_idx - 1;
                linear_idx = sub2ind(size(dd), row, global_col_idx);
                estimated_non_zero_indices = [estimated_non_zero_indices; linear_idx];
                
                im_bits_est = de2bi(detected_local_idx - 1, bits_per_IM_data, 'left-msb');
                qam_bits_est = qamdemod(detected_symbol, cen_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
                estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
            end
        end

        % --- B: 处理保护间隔区域 (使用相同的联合ML检测逻辑) ---
        for row = guard_rows
            i_col = 1;
            while i_col <= num_guard_cols
                current_group_size = min(g, num_guard_cols - i_col + 1);
                if current_group_size  > 1
                    global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
                    group_est = dd_est_matrix(row, global_col_indices) / scale_guard;
                    
                    max_metric = -inf;
                    best_local_idx = -1;
                    best_symbol = 0;

                    % 遍历所有可能激活的位置
                    for k = 1:current_group_size
                        % 遍历所有可能的QAM星座点
                        for q = 1:guard_modu
                            current_symbol = guard_constellation(q);
                            metric = -abs(group_est(k) - current_symbol)^2;
                            for m = 1:current_group_size
                                if m ~= k
                                    metric = metric - abs(group_est(m))^2;
                                end
                            end
                            if metric > max_metric
                                max_metric = metric;
                                best_local_idx = k;
                                best_symbol = current_symbol;
                            end
                        end
                    end

                    detected_local_idx = best_local_idx;
                    detected_symbol = best_symbol;
                    
                    global_col_idx = global_col_indices(detected_local_idx);
                    linear_idx = sub2ind(size(dd), row, global_col_idx);
                    estimated_non_zero_indices = [estimated_non_zero_indices; linear_idx];
                    
                    % 解码比特流 (注意处理组大小不是g的情况)
                    if current_group_size == g
                        im_bits_est = de2bi(detected_local_idx - 1, bits_per_IM_data, 'left-msb');
                    else % 独热编码
                        im_bits_est = zeros(1, current_group_size);
                        im_bits_est(detected_local_idx) = 1;
                    end
                    qam_bits_est = qamdemod(detected_symbol, guard_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
                    estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
                end
                i_col = i_col + current_group_size;
            end
        end

        %
        % [rows, cols] = find(dd);
        % sorted_positions0 = sortrows([rows, cols], 1);
        % 
        % estimated_non_zero_indices = sort(estimated_non_zero_indices);
        % [estimated_rows, estimated_cols] = ind2sub([M, N], estimated_non_zero_indices);
        % estimated_positions_matrix = [estimated_rows, estimated_cols];
        % sorted_positions1 = sortrows(estimated_positions_matrix, 1);
        % missed_positions = setdiff(sorted_positions0, sorted_positions1, 'rows');

        [num_errors_this_frame, ~] = biterr(original_bits_stream, estimated_bits_stream);
        num_bits_this_frame = length(original_bits_stream);
        frame_errors(f) = num_errors_this_frame;
        frame_bits(f) = num_bits_this_frame;
    end
    total_bit_errors = sum(frame_errors);
    total_bits_transmitted = sum(frame_bits);
    ber(i_snr) =total_bit_errors / total_bits_transmitted;% 计算当前信噪比下的平均BER
end

figure();
semilogy(SNR, ber, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('信噪比 (dB)');ylabel('比特误码率 (BER)');title('MMSE检测后系统BER性能曲线');legend('IM-proposed');

function output = zeta_N(k, N)
tolerance = 1e-10; % 可调容差
if (k == N||k==0)
    output = 1;
else
    exp_term = exp(-1i * pi * (N - 1) * k / N);
    numerator = sin(pi * k); % 分子 sin(pi * k)
    if abs(numerator) < tolerance
        numerator = 0;
    end
    denominator = sin(pi * k / N); % 分母 sin(pi * k / N)
    output = exp_term * (numerator / (N * denominator));
end
end
