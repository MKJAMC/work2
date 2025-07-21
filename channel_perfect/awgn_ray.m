% % % %% BPSK over Rayleigh Fading Channel with Perfect CSI
% % % %  ----------------------------------------------------
% % % %  此代码用于展示在瑞利衰落信道下，即使拥有完美的信道信息，
% % % %  其性能也远差于AWGN信道，并以此来对比您之前代码中的理想化结果。
% % % 
% % % clc;
% % % clear;
% % % close all;
% % % 
% % % %% 1. 参数定义
% % % % =================================
% % % N_bits = 1e6;       % 仿真总比特数 (为了在低BER时获得平滑曲线，数量要足够大)
% % % SNR_dB_range = 0:2:20; % 信噪比范围 (dB)
% % % M = 2;              % 调制阶数 (BPSK)
% % % 
% % % % 初始化结果存储向量
% % % ber_simulated = zeros(1, length(SNR_dB_range));
% % % ber_theoretical_rayleigh = zeros(1, length(SNR_dB_range));
% % % ber_theoretical_awgn = zeros(1, length(SNR_dB_range));
% % % 
% % % 
% % % %% 2. 仿真主循环
% % % % =================================
% % % for idx = 1:length(SNR_dB_range)
% % % 
% % %     snr_db = SNR_dB_range(idx);
% % %     fprintf('正在仿真 SNR = %d dB...\n', snr_db);
% % % 
% % %     % --- 数据生成与调制 ---
% % %     bits = randi([0 1], 1, N_bits);         % 生成随机比特流
% % %     s = 2*bits - 1;                         % BPSK调制 (0 -> -1, 1 -> +1)
% % % 
% % %     % --- 信道建模 ---
% % %     % 1. 生成瑞利衰落信道系数 h
% % %     % h 的模服从瑞利分布, 功率为1
% % %     h = (randn(1, N_bits) + 1i*randn(1, N_bits)) / sqrt(2);
% % % 
% % %     % 2. 信号通过衰落信道
% % %     s_faded = h .* s;
% % % 
% % %     % 3. 添加高斯白噪声 (AWGN)
% % %     % 将dB信噪比转换为线性值 (Eb/N0)
% % %     % 对于BPSK, 符号能量Es = 比特能量Eb = 1
% % %     snr_linear = 10^(snr_db/10);
% % % 
% % %     % 计算噪声功率
% % %     noise_power = 1 / snr_linear;
% % % 
% % %     % 生成复数噪声
% % %     noise = sqrt(noise_power/2) * (randn(1, N_bits) + 1i*randn(1, N_bits));
% % % 
% % %     % 接收信号
% % %     y = s_faded + noise;
% % % 
% % %     % --- 接收机处理 ---
% % %     % 1. 信道均衡 (核心步骤)
% % %     %  *** 假设拥有完美的信道信息h ***
% % %     y_eq = y ./ h;
% % % 
% % %     % 2. BPSK解调
% % %     % 直接对均衡后信号的实部进行判决
% % %     received_bits = real(y_eq) > 0;
% % % 
% % %     % --- 性能计算 ---
% % %     % 计算误码数
% % %     n_errors = sum(received_bits ~= bits);
% % % 
% % %     % 计算并存储当前信噪比下的误码率
% % %     ber_simulated(idx) = n_errors / N_bits;
% % % 
% % %     % --- 计算理论值用于对比 ---
% % %     % 理论AWGN误码率
% % %     ber_theoretical_awgn(idx) = qfunc(sqrt(2*snr_linear));
% % % 
% % %     % 理论瑞利衰落误码率 (完美CSI下)
% % %     ber_theoretical_rayleigh(idx) = 0.5 * (1 - sqrt(snr_linear / (1 + snr_linear)));
% % % 
% % % end
% % % 
% % % 
% % % %% 3. 结果绘图
% % % % =================================
% % % figure;
% % % semilogy(SNR_dB_range, ber_theoretical_awgn, 'b-', 'LineWidth', 2);
% % % hold on;
% % % semilogy(SNR_dB_range, ber_theoretical_rayleigh, 'r--', 'LineWidth', 2);
% % % semilogy(SNR_dB_range, ber_simulated, 'ro', 'MarkerSize', 8);
% % % grid on;
% % % xlabel('信噪比 (SNR / dB)');
% % % ylabel('比特误码率 (BER)');
% % % title('BPSK在瑞利衰落信道与AWGN信道下的性能对比');
% % % legend('理论AWGN', '理论瑞利 (完美CSI)', '仿真瑞利 (完美CSI)');
% % % axis([min(SNR_dB_range) max(SNR_dB_range) 1e-5 1]);
% % 
% % 
% % %% BPSK over Rayleigh Fading Channel with Perfect CSI (Fixed Noise Power)
% % %  --------------------------------------------------------------------
% % %  此代码遵循另一种仿真逻辑：固定噪声功率为1，通过调整信号功率
% % %  来满足给定的信噪比(SNR)要求。其最终结果与前一版代码完全相同。
% % 
% % clc;
% % clear;
% % close all;
% % 
% % %% 1. 参数定义
% % % =================================
% % N_bits = 1e6;       % 仿真总比特数
% % SNR_dB_range = 0:2:20; % 信噪比范围 (dB)
% % M = 2;              % 调制阶数 (BPSK)
% % 
% % % 初始化结果存储向量
% % ber_simulated = zeros(1, length(SNR_dB_range));
% % ber_theoretical_rayleigh = zeros(1, length(SNR_dB_range));
% % ber_theoretical_awgn = zeros(1, length(SNR_dB_range));
% % 
% % 
% % %% 2. 仿真主循环
% % % =================================
% % for idx = 1:length(SNR_dB_range)
% % 
% %     snr_db = SNR_dB_range(idx);
% %     fprintf('正在仿真 SNR = %d dB...\n', snr_db);
% % 
% %     % 将dB信噪比转换为线性值 (Eb/N0)
% %     snr_linear = 10^(snr_db/10);
% % 
% %     % --- 数据生成与调制 (核心改动部分) ---
% %     % 1. 固定噪声总功率 N0 为 1
% %     N0 = 1;
% % 
% %     % 2. 根据 SNR = Eb/N0，计算所需的比特能量 Eb
% %     Eb = snr_linear * N0; % 因为 N0 = 1, 所以 Eb = snr_linear
% % 
% %     % 3. 计算BPSK符号所需的幅度 A (因为 Eb = A^2)
% %     A = sqrt(Eb);
% % 
% %     % 4. 生成比特并进行幅度调整后的BPSK调制
% %     bits = randi([0 1], 1, N_bits);
% %     s = A * (2*bits - 1);  % 符号为 -A 和 +A
% % 
% %     % --- 信道建模 ---
% %     % 1. 生成瑞利衰落信道系数 h
% %     h = (randn(1, N_bits) + 1i*randn(1, N_bits)) / sqrt(2);
% % 
% %     % 2. 信号通过衰落信道
% %     s_faded = h .* s;
% % 
% %     % 3. 添加高斯白噪声 (核心改动部分)
% %     % 噪声总功率固定为1
% %     noise_power = N0; % 等于 1
% %     noise = sqrt(noise_power/2) * (randn(1, N_bits) + 1i*randn(1, N_bits));
% % 
% %     % 接收信号
% %     y = s_faded + noise;
% % 
% %     % --- 接收机处理 ---
% %     % (这部分逻辑与之前完全一样)
% %     % 1. 信道均衡 (假设拥有完美的信道信息h)
% %     y_eq = y ./ h;
% % 
% %     % 2. BPSK解调
% %     received_bits = real(y_eq) > 0;
% % 
% %     % --- 性能计算 ---
% %     n_errors = sum(received_bits ~= bits);
% %     ber_simulated(idx) = n_errors / N_bits;
% % 
% %     % --- 计算理论值用于对比 ---
% %     % (理论公式不变，因为它只依赖于Eb/N0这个比值)
% %     ber_theoretical_awgn(idx) = qfunc(sqrt(2*snr_linear));
% %     ber_theoretical_rayleigh(idx) = 0.5 * (1 - sqrt(snr_linear / (1 + snr_linear)));
% % 
% % end
% % 
% % 
% % %% 3. 结果绘图
% % % =================================
% % figure;
% % semilogy(SNR_dB_range, ber_theoretical_awgn, 'b-', 'LineWidth', 2);
% % hold on;
% % semilogy(SNR_dB_range, ber_theoretical_rayleigh, 'r--', 'LineWidth', 2);
% % semilogy(SNR_dB_range, ber_simulated, 'ro', 'MarkerSize', 8);
% % grid on;
% % xlabel('信噪比 (SNR / dB)');
% % ylabel('比特误码率 (BER)');
% % title('BPSK在瑞利信道下的性能对比 (固定噪声功率)');
% % legend('理论AWGN', '理论瑞利 (完美CSI)', '仿真瑞利 (完美CSI)');
% % axis([min(SNR_dB_range) max(SNR_dB_range) 1e-5 1]);
% 
% 
% %% QPSK over Rayleigh Fading Channel with Perfect CSI (Fixed Noise Power)
% %  --------------------------------------------------------------------
% %  此代码将之前的BPSK仿真修改为QPSK。
% %  关键变化在于处理每个符号2比特的能量关系以及调制解调过程。
% 
% clc;
% clear;
% close all;
% 
% %% 1. 参数定义
% % =================================
% N_bits = 2e6;       % 仿真总比特数 (QPSK每个符号2比特，建议为偶数)
% M = 4;              % 调制阶数 (QPSK)
% k = log2(M);        % 每个符号携带的比特数 (k=2)
% 
% SNR_dB_range = 0:2:20; % 信噪比范围 (dB)
% 
% % 确保比特总数为偶数
% N_bits = ceil(N_bits/k)*k;
% 
% % 初始化结果存储向量
% ber_simulated = zeros(1, length(SNR_dB_range));
% ber_theoretical_rayleigh = zeros(1, length(SNR_dB_range));
% ber_theoretical_awgn = zeros(1, length(SNR_dB_range));
% 
% 
% %% 2. 仿真主循环
% % =================================
% for idx = 1:length(SNR_dB_range)
% 
%     snr_db = SNR_dB_range(idx);
%     fprintf('正在仿真 SNR = %d dB...\n', snr_db);
% 
%     % 将dB信噪比转换为线性值 (Eb/N0)
%     snr_linear = 10^(snr_db/10);
% 
%     % --- 数据生成与调制 (核心改动部分) ---
%     % 1. 固定噪声总功率 N0 为 1
%     N0 = 1;
% 
%     % 2. 根据 SNR = Eb/N0，计算所需的比特能量 Eb
%     Eb = snr_linear * N0;
% 
%     % 3. 计算QPSK符号所需的符号能量 Es (Es = 2 * Eb)
%     Es = k * Eb;
% 
%     % 4. 计算信号幅度 A (A^2 = Es)
%     A = sqrt(Es);
% 
%     % 5. 生成随机比特流
%     bits = randi([0 1], N_bits, 1);
% 
%     % 6. 使用qammod进行QPSK调制
%     %    'UnitAveragePower', true 使得输出符号的平均能量为1
%     %    它会自动进行格雷码映射
%     s_base = qammod(bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
% 
%     % 7. 按计算出的幅度A对信号进行能量缩放
%     s = A * s_base;
% 
%     % --- 信道建模 ---
%     N_symbols = N_bits / k;
%     % 1. 生成瑞利衰落信道系数 h (每个符号一个)
%     h = (randn(N_symbols, 1) + 1i*randn(N_symbols, 1)) / sqrt(2);
% 
%     % 2. 信号通过衰落信道
%     s_faded = h .* s;
% 
%     % 3. 添加高斯白噪声
%     noise_power = N0; % 等于 1
%     noise = sqrt(noise_power/2) * (randn(N_symbols, 1) + 1i*randn(N_symbols, 1));
% 
%     % 接收信号
%     y = s_faded + noise;
% 
%     % --- 接收机处理 ---
%     % 1. 信道均衡 (完美CSI)
%     y_eq = y ./ h;
% 
%     % 2. QPSK解调 (核心改动部分)
%     %    解调前，需要将接收符号的能量归一化回1，即除以幅度A
%     %    这样才能匹配'UnitAveragePower'为true的解调器
%     received_bits = qamdemod(y_eq / A, M, 'OutputType', 'bit', 'UnitAveragePower', true);
% 
%     % --- 性能计算 ---
%     n_errors = sum(received_bits ~= bits);
%     ber_simulated(idx) = n_errors / N_bits;
% 
%     % --- 计算理论值用于对比 ---
%     % QPSK的理论误码率(格雷码)近似于BPSK的公式
%     ber_theoretical_awgn(idx) = qfunc(sqrt(2*snr_linear));
%     ber_theoretical_rayleigh(idx) = 0.5 * (1 - sqrt(snr_linear / (1 + snr_linear)));
% 
% end
% 
% 
% %% 3. 结果绘图
% % =================================
% figure;
% semilogy(SNR_dB_range, ber_theoretical_awgn, 'b-', 'LineWidth', 2);
% hold on;
% semilogy(SNR_dB_range, ber_theoretical_rayleigh, 'r--', 'LineWidth', 2);
% semilogy(SNR_dB_range, ber_simulated, 'ms', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
% grid on;
% xlabel('每比特信噪比 (E_b/N_0 / dB)');
% ylabel('比特误码率 (BER)');
% title('QPSK在瑞利衰落信道与AWGN信道下的性能对比');
% legend('理论AWGN', '理论瑞利 (近似)', '仿真QPSK (瑞利)');
% axis([min(SNR_dB_range) max(SNR_dB_range) 1e-5 1]);

% 1. EVA信道模型的相对功率 (dB) 
relative_power_dB = [0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];

% 步骤一: 将dB转换为线性值
relative_power_linear = 10.^(relative_power_dB / 10);

% 步骤二: 功率归一化
total_power = sum(relative_power_linear);
normalized_variances = relative_power_linear / total_power;%计算每个系数的百分比

% 检查：所有方差之和应约等于1
disp(['总归一化功率: ', num2str(sum(normalized_variances))]);

% `normalized_variances` 这个数组现在包含了每条路径的方差 sigma_i^2
% 可以在仿真循环中用来生成 h_i
disp('每条路径的归一化方差 (sigma_i^2):');
disp(normalized_variances);

% 步骤三: 在仿真循环中生成随机信道增益的示例
num_paths = length(normalized_variances);
h = zeros(1, num_paths); % 存储当前时刻的信道增益

for i = 1:num_paths
    sigma2_i = normalized_variances(i);
    % 生成一个随机的信道增益
    h(i) = sqrt(sigma2_i / 2) * (randn() + 1j * randn());
end

disp('一次仿真迭代生成的随机信道增益 h_i:');
disp(h);