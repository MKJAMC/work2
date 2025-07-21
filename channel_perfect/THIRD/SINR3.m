clc;
clear;
close all;

%% 1. 参数定义 (完全遵循方案三的设定)
M = 64;
N = 32;
l_max = 20;
k_max = 2;
g = 4; % IM分组参数

% 信道参数
fc=64e9; delta_f=120e3; u_max=120*1000/3600; c=3e8;
delay=[30,150]*1e-9;
delay_tap0=delay*M*delta_f;
li=round(delay_tap0);
doppler_max=fc*u_max/c;
theta0 = -pi + 2*pi*rand(1, length(delay));
doppler=doppler_max*cos(theta0);
ki=doppler*N*(1/delta_f);
h_p_db=[1.5,1.4];
relative_power_linear = 10.^((-1) * h_p_db/10);
total_power = sum(relative_power_linear);
h_p = relative_power_linear / total_power;
h_exp = exp( 1j*2 * pi * rand(1, length(h_p_db)));
P = length(delay);

% 方案三特有参数
cen_modu = 4;
guard_modu = 2;
bits_per_qam_symbol = log2(cen_modu);
bits_per_guard_symbol=log2(guard_modu);
%% 
guard_data_power_factor = 1;

% 分析的SNR点
SNR_point = 10; % dB
base_power_persym = 10.^(SNR_point/10);

%% 2. 功率归一化 (以方案二为基准)
% 计算方案二和方案三的激活符号数
ans2 = (M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
mo = mod((N-4*k_max-1), g);
q1 = floor((N-4*k_max-1) / g);
cen_num = (N/g)*(M-(1+2*l_max));
bit_row_num = (1+q1)*(2*l_max+1);
ans3 = cen_num + bit_row_num;

% 计算方案三的每符号平均功率
total_frame_power = base_power_persym * ans2;
factor_s3 = total_frame_power / ans3;

fprintf('================== 功率设定 ==================\n');
fprintf('总目标功率 (基于方案二) = %.4f\n', total_frame_power);
fprintf('方案三激活符号数 = %d\n', ans3);
fprintf('方案三每符号平均功率 (factor) = %.4f\n', factor_s3);
fprintf('==============================================\n');

%% 3. 生成单次快照的方案三 (OTFS-IM) 信号
disp('正在生成方案三 (IM) 的发射信号...');
% --- 功率分配计算 ---
num_central_symbols = cen_num;
num_guard_symbols = ans3 - num_central_symbols;
denominator = num_central_symbols + num_guard_symbols * guard_data_power_factor;
power_per_central_symbol = total_frame_power / denominator;
power_per_guard_symbol = power_per_central_symbol * guard_data_power_factor;
scale_central = sqrt(power_per_central_symbol);
scale_guard = sqrt(power_per_guard_symbol);

% --- 生成IM信号 (仅生成一帧用于分析) ---
dd = zeros(M, N);
% *** MODIFICATION START: Initialize index trackers ***
central_indices_vec = [];
guard_indices_vec = [];
% *** MODIFICATION END ***

data_rows = (l_max + 2):(M - l_max);
num_groups_data = N / g;
bits_per_IM_data = log2(nchoosek(g, 1));
% A: 处理中央数据区域
for row = data_rows
    for j = 1:num_groups_data
        im_bits = randi([0 1], 1, bits_per_IM_data);
        qam_bits = randi([0 1], 1, bits_per_qam_symbol);
        active_local_idx = bi2de(im_bits, 'left-msb') + 1;
        symbol_int = bi2de(qam_bits, 'left-msb');
        qam_symbol = qammod(symbol_int, cen_modu, 'UnitAveragePower', true);
        start_col = (j - 1) * g + 1;
        active_col = start_col + active_local_idx - 1;
        dd(row, active_col) = qam_symbol;
        
        % *** MODIFICATION START: Track central symbol index ***
        linear_idx = (active_col - 1) * M + row;
        central_indices_vec = [central_indices_vec; linear_idx];
        % *** MODIFICATION END ***
    end
end
% B: 处理保护间隔区域
guard_rows = [1:(l_max + 1), (M - l_max + 1):M];
guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
num_guard_cols = length(guard_cols_range);
for row = guard_rows
    i_col = 1;
    while i_col <= num_guard_cols
        current_group_size = min(g, num_guard_cols - i_col + 1);
        if current_group_size > 1
            if current_group_size == g
                im_bits = randi([0 1], 1, bits_per_IM_data);
                active_local_idx = bi2de(im_bits, 'left-msb') + 1;
            else
                active_local_idx = randi(current_group_size);
            end
            qam_bits = randi([0 1], 1, bits_per_guard_symbol);
            symbol_int = bi2de(qam_bits, 'left-msb');
            qam_symbol = qammod(symbol_int, guard_modu, 'UnitAveragePower', true);
            global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
            active_col = global_col_indices(active_local_idx);
            dd(row, active_col) = qam_symbol;
            
            % *** MODIFICATION START: Track guard symbol index ***
            linear_idx = (active_col - 1) * M + row;
            guard_indices_vec = [guard_indices_vec; linear_idx];
            % *** MODIFICATION END ***
        end
        i_col = i_col + current_group_size;
    end
end
% 应用功率缩放因子
dd(data_rows,:) = dd(data_rows,:) * scale_central;
guard_data_indices_rows = [1:(l_max+1), (M-l_max+1):M];
dd(guard_data_indices_rows,:) = dd(guard_data_indices_rows,:) * scale_guard;
dd_vec = dd(:); % 最终的发射向量
disp('发射信号生成完毕。');

%% 4. 核心计算 (构建H_eff, G, GH)
disp('正在执行核心计算...');
hw=zeros(M,N);
for l=0:M-1, for k=0:N-1, for i=1:P, theta=exp(1j*2*pi*ki(i)*(li(i))/(M*N)); delta_term = (l == li(i)); hw(l+1,k+1)=hw(l+1,k+1)+h_p(i)*h_exp(i)*delta_term*zeta_N(k-ki(i),N).*theta; end, end, end
H_eff = zeros(M*N, M*N);
for l=0:M-1, for k=0:N-1, row_idx=k*M+(l+1); for l_prime=0:M-1, for k_prime=0:N-1, col_idx=k_prime*M+(l_prime+1); h_l=mod(l-l_prime,M); h_k=mod(k-k_prime,N); H_eff(row_idx,col_idx)=hw(h_l+1,h_k+1); end, end, end, end
heff = H_eff'; Ie = eye(size(H_eff,2));

% 使用方案三的平均每符号功率 factor_s3 来构建其MMSE均衡器
G = (heff * ((H_eff * heff + (1/factor_s3) * Ie))^(-1));

GH = G * H_eff;
GG_h = G * G';
disp('核心矩阵计算完成。');

%% 5. 计算瞬时S, I, N以及SINR (采用第一段代码的公式)
disp('正在计算瞬时SINR...');
total_received_vec = GH * dd_vec;
signal_only_vec = diag(diag(GH)) * dd_vec;
interference_vec = total_received_vec - signal_only_vec;

signal_power_inst = abs(signal_only_vec).^2;
interference_power_inst = abs(interference_vec).^2;
noise_power_inst = diag(GG_h); % 均衡后的噪声功率

% 瞬时SINR (只对非零发射符号有意义)
sinr_inst = signal_power_inst ./ (interference_power_inst + noise_power_inst + 1e-20);
disp('SINR计算完成。');

%% 6. 结果分析
% *** MODIFICATION START: Separate SINR analysis ***

% A: 分析中央区域符号
% 提取中央区域激活符号对应的瞬时SINR
sinr_central = sinr_inst(central_indices_vec);
% 计算平均值
avg_sinr_central_linear = mean(sinr_central);
avg_sinr_central_db = 10*log10(avg_sinr_central_linear);

% B: 分析保护区域符号
% 提取保护区域激活符号对应的瞬时SINR
sinr_guard = sinr_inst(guard_indices_vec);
% 计算平均值
avg_sinr_guard_linear = mean(sinr_guard);
avg_sinr_guard_db = 10*log10(avg_sinr_guard_linear);


fprintf('\n================== 方案三结果分析 ==================\n');
fprintf('中央区域: 激活符号数 = %d\n', length(central_indices_vec));
fprintf('中央区域: 平均SINR = %.4f (线性值)\n', avg_sinr_central_linear);
fprintf('中央区域: 平均SINR = %.4f (dB)\n\n', avg_sinr_central_db);

fprintf('保护区域: 激活符号数 = %d\n', length(guard_indices_vec));
fprintf('保护区域: 平均SINR = %.4f (线性值)\n', avg_sinr_guard_linear);
fprintf('保护区域: 平均SINR = %.4f (dB)\n', avg_sinr_guard_db);
fprintf('====================================================\n');

% *** MODIFICATION END ***

%% 7. zeta_N 辅助函数
function output = zeta_N(k, N)
    tolerance = 1e-10;
    if (abs(k) < tolerance || abs(k-N) < tolerance), output = 1;
    else, exp_term = exp(-1i * pi * (N-1)*k/N); numerator = sin(pi*k);
        if abs(numerator)<tolerance, numerator=0; end; denominator = sin(pi*k/N);
        if abs(denominator)<tolerance, output=0; else, output=exp_term*(numerator/(N*denominator)); end;
    end
end