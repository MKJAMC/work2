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
theta0= -pi + 2*pi*rand(1, length(delay));
doppler=doppler_max*cos(theta0) ;
ki=doppler*N*(1/delta_f);
k_max0=doppler_max*N*(1/delta_f);
% EVA信道，时延为0的路径是增益最大的，以后依次递减
h_p_db=[1.5,1.4,3.6,0.6,9.1,7,12,16.9];
relative_power_linear = 10.^((-1) * h_p_db/10);%db转换为功率
total_power = sum(relative_power_linear);
h_p = relative_power_linear / total_power;%每一个增益所占据的百分比
h_exp=exp( 1j*2 * pi * rand(1, length(h_p_db)));

l_max=20;k_max=2;
P=length(delay);%路径数
cen_modu = 2;
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
bit_row_num=(1+q1)*(2*l_max+1);%guard个数
cen_num=(N/g)*(M-(1+2*l_max));
ans3=cen_num+bit_row_num;
bit=(log2(cen_modu)+log2(nchoosek(g,1)))*(cen_num)+(2*l_max+1)*bit_row;

%% 符号检测SNR
SNR=10:2:16;
iter=[4e4,4e4,4e5,4e5];
% SNR=18:2:20;iter=[4e5,4e6];
power_persym=10.^(SNR/10);%每个符号的功率
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
power=power_persym*ans2;%ans2为基准，求得总功率，方案二的功率为1
ber = zeros(1, length(SNR));
factor=power/ans3;

for i_snr=1:length(SNR)
    ber_iter = zeros(1, iter(i_snr));%一个frame错误bit数
    num_frames =ceil(iter(i_snr)/bit);
    frame_im_errors = zeros(1, num_frames);
    frame_qam_errors = zeros(1, num_frames);
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
    % 3.4 计算幅度缩放因子 (幅度是功率的平方根
    scale_central = sqrt(power_per_central_symbol);
    scale_guard   = sqrt(power_per_guard_symbol);

    for f=1:num_frames
        disp(f)
        original_im_bits = [];        original_qam_bits = [];
        %%  发射端：基于信息比特生成OTFS-IM信号
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

                original_im_bits = [original_im_bits, im_bits];                original_qam_bits = [original_qam_bits, qam_bits];

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
            while i_col <= num_guard_cols     % 计算当前组的大小，最大不超过 g
                current_group_size = min(g, num_guard_cols - i_col + 1);
                if current_group_size > 1  %最后余数为1，则不存在索引调制
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
                    original_im_bits = [original_im_bits, im_bits];
                    original_qam_bits = [original_qam_bits, qam_bits];

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
        dd(guard_data_indices,:)   = dd(guard_data_indices,:)* scale_guard;


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

    

        %% H_eff等效信道
        y_vec = reshape(Y, M*N, 1);
        dd_vec = reshape(dd, M*N, 1);
        H_eff = zeros(M*N, M*N);
        for l = 0:M-1       % 遍历矩阵的行
            for k = 0:N-1   % 遍历矩阵的列
                % 计算 Y(l+1, k+1) 在列主序向量 y_vec 中的索引
                row_idx = k*M + (l+1);
                for l_prime = 0:M-1
                    for k_prime = 0:N-1
                        % 计算 dd(l_prime+1, k_prime+1) 在列主序向量 dd_vec 中的索引
                        col_idx = k_prime*M + (l_prime+1);
                        h_l = mod(l - l_prime, M);
                        h_k = mod(k - k_prime, N);
                        H_eff(row_idx, col_idx) = hw(h_l + 1, h_k + 1);
                    end
                end
            end
        end
        %% 添加噪声
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声
        y=H_eff*dd_vec+noise;
       
      
        % MMSE 检测 
        heff=H_eff';
        Ie=eye(size(H_eff, 2));
        H_lmmse=(heff*((H_eff*heff+(1/factor(i_snr))*Ie))^(-1));%平均snr每个符号的
        dd_est_matrix=H_lmmse*y;
        dd_est_matrix= reshape(dd_est_matrix, M, []);

        estimated_bits_stream=[];
        estimated_im_bits=[];estimated_qam_bits=[];

        % ... 此前的代码已经计算出了估计矩阵 dd_est_matrix ...

        % 初始化一个空数组，用于存储在接收端检测出的比特流
        estimated_bits_stream = [];

        % % --- A: 对中心数据区域进行检测 ---
        % % 我们遍历发射端已知的网格结构，来检测传输的比特。
        % for row = data_rows
        %     for j = 1:num_groups_data
        %         % 定义当前索引调制(IM)组所占的列
        %         start_col = (j - 1) * g + 1;
        %         group_cols = start_col : start_col + g - 1;
        %         % 从原始发送矩阵 dd 中提取当前组的符号，以找到真实激活位置
        %         original_group_symbols = dd(row, group_cols);
        %         %% 【关键修改】: 直接使用真实位置，不再搜索最大值
        %         [~, true_local_idx] = max(abs(original_group_symbols));
        %         im_bits_est = de2bi(true_local_idx - 1, bits_per_IM_data, 'left-msb');
        %         group_est_symbols = dd_est_matrix(row, group_cols) / scale_central;
        %         symbol_to_demod = group_est_symbols(true_local_idx);
        %         qam_bits_est = qamdemod(symbol_to_demod, cen_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
        %         estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
        %     end
        % end
        % 
        % % --- B: 对保护间隔区域进行检测 ---
        % for row = guard_rows
        %     i_col = 1;
        %     while i_col <= num_guard_cols
        %         current_group_size = min(g, num_guard_cols - i_col + 1);
        %         if current_group_size > 1
        %             global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
        %             %% 【关键修改】: 直接从原始 dd 矩阵中找到真实激活位置
        %             original_group_symbols = dd(row, global_col_indices);
        %             [~, true_local_idx] = max(abs(original_group_symbols));
        %             % 1. 索引检测 (理想化)
        %             if current_group_size == g
        %                 num_im_bits = floor(log2(nchoosek(g, 1)));
        %                 im_bits_est = de2bi(true_local_idx - 1, num_im_bits, 'left-msb');
        %             else
        %                 im_bits_est = zeros(1, current_group_size);
        %                 im_bits_est(true_local_idx) = 1;
        %             end
        %             % 2. 符号检测
        %             group_est_symbols = dd_est_matrix(row, global_col_indices) / scale_guard;
        %             symbol_to_demod = group_est_symbols(true_local_idx);
        %             qam_bits_est = qamdemod(symbol_to_demod, guard_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
        %             estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
        %         end
        %         i_col = i_col + current_group_size;
        %     end
        % end

        % --- A: 处理中央数据区域 ---
        % 检测中央数据区域
        for row = data_rows
            for j = 1:num_groups_data
                start_col = (j - 1) * g + 1;
                group_est = dd_est_matrix(row, start_col : start_col + g - 1) / scale_central;
                [~, detected_local_idx] = max(abs(group_est));
                detected_symbol = group_est(detected_local_idx);
                im_bits_est = de2bi(detected_local_idx - 1, bits_per_IM_data, 'left-msb');
                qam_bits_est = qamdemod(detected_symbol, cen_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
                estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
                %
                estimated_im_bits = [estimated_im_bits, im_bits_est];
                estimated_qam_bits = [estimated_qam_bits, qam_bits_est(:).'];
            end
        end
        len1=length(estimated_bits_stream);
        [~,err1]=biterr(original_bits_stream(1:len1), estimated_bits_stream);

        % --- B: 处理保护间隔区域 ---
        for row = guard_rows
            i_col = 1;
            while i_col <= num_guard_cols
                current_group_size = min(g, num_guard_cols - i_col + 1);
                if current_group_size  > 1 %与发射端设计呼应
                    if current_group_size == g
                        global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
                        group_est = dd_est_matrix(row, global_col_indices) / scale_guard;
                        [~, detected_local_idx] = max(abs(group_est));
                        detected_symbol = group_est(detected_local_idx);
                        im_bits_est = de2bi(detected_local_idx - 1, bits_per_IM_data, 'left-msb');
                        qam_bits_est = qamdemod(detected_symbol, guard_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
                        estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
                        %
                        estimated_im_bits = [estimated_im_bits, im_bits_est];
                        estimated_qam_bits = [estimated_qam_bits, qam_bits_est(:).'];
                    else
                        global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
                        group_est = dd_est_matrix(row, global_col_indices) / scale_guard;
                        [~, detected_local_idx] = max(abs(group_est));
                        detected_symbol = group_est(detected_local_idx);
                        im_bits_est = zeros(1, current_group_size);
                        im_bits_est(detected_local_idx) = 1;
                        qam_bits_est = qamdemod(detected_symbol, guard_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
                        estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
                        %
                        estimated_im_bits = [estimated_im_bits, im_bits_est];
                        estimated_qam_bits = [estimated_qam_bits, qam_bits_est(:).'];
                    end
                end
                i_col = i_col + current_group_size;
            end
        end
        [~,err2]=biterr(original_bits_stream(len1+1:length(estimated_bits_stream)), estimated_bits_stream(len1+1:length(estimated_bits_stream)));
        %% 调试
        % A=reshape(original_bits_stream(len1+1:length(estimated_bits_stream)),bit_row,[]);
        % B=reshape(estimated_bits_stream(len1+1:length(estimated_bits_stream)),bit_row,[]);
        % IM1=q1*(log2(guard_modu)+log2(nchoosek(g,1)));
        % 
        % [num1,~]=biterr(A(1:IM1,:),B(1:IM1,:));%正常调制
        % [num2,~]=biterr(A(IM1+1:bit_row,:),B(IM1+1:bit_row,:));%独热编码
      
        %
        [im_err_count, ~] = biterr(original_im_bits, estimated_im_bits);
        [qam_err_count, ~] = biterr(original_qam_bits, estimated_qam_bits);
        frame_im_errors(f) = im_err_count;
        frame_qam_errors(f) = qam_err_count;
        frame_total_im_bits(f) = length(original_im_bits);
        frame_total_qam_bits(f) = length(original_qam_bits);

        [num_errors_this_frame, ~] = biterr(original_bits_stream, estimated_bits_stream);
        num_bits_this_frame = length(original_bits_stream);
        frame_errors(f) = num_errors_this_frame;
        frame_bits(f) = num_bits_this_frame;

    end
    total_im_errors = sum(frame_im_errors);    total_qam_errors = sum(frame_qam_errors);
    total_im_bits = sum(frame_total_im_bits);    total_qam_bits = sum(frame_total_qam_bits);
    ber_im(i_snr) = total_im_errors / total_im_bits; ber_qam(i_snr) = total_qam_errors / total_qam_bits;

    total_bit_errors = sum(frame_errors);
    total_bits_transmitted = sum(frame_bits);
    ber(i_snr) =total_bit_errors / total_bits_transmitted;% 计算当前信噪比下的平均BER
end


figure();
semilogy(SNR, ber, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('信噪比 (dB)');ylabel('比特误码率 (BER)');title('MMSE检测后系统BER性能曲线');legend('MMSE-BER');

function output = zeta_N(k, N)
%这个函数属于输入一个k，得到一个数值
tolerance = 1e-10;  % 可调容差
if (k == N||k==0)
    output = 1;
else
    exp_term = exp(-1i * pi * (N - 1) * k / N);
    numerator = sin(pi * k);  % 分子 sin(pi * k)
    if abs(numerator) < tolerance
        numerator = 0;
    end
    denominator = sin(pi * k / N);  % 分母 sin(pi * k / N)
    output = exp_term * (numerator / (N * denominator));
end
end

