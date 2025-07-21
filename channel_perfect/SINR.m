clc;
clear;
close all;

%% 1. 参数定义
M = 64; 
N = 32; 
l_max = 20; 
k_max = 2;

% 信道参数 (简化为2条路径)
fc=64e9; delta_f=120e3; u_max=120*1000/3600; c=3e8;
delay=[30,150]*1e-9;
delay_tap0=delay*M*delta_f; 
li=round(delay_tap0); % li 将会是 [0, 1]
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

% 分析的SNR点
SNR_point = 10; % dB
base_power_persym = 10.^(SNR_point/10); 

% --- 修正与明确功率归一化逻辑 ---
% !!! 目标：以方案二的总功率为基准 !!!
ans1 = (M-2*l_max-1)*N;
ans2 = (M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
ans3=246+184;
% 1. 方案二的单符号功率即为基准功率
factor_s2 = base_power_persym;
% 2. 计算出帧的总功率
total_frame_power = factor_s2 * ans2;
% 3. 计算方案一为了达到相同总功率，其单符号应有的功率
factor_s1 = total_frame_power / ans1;
factor_s3=total_frame_power/ans3;
% fprintf('功率设定: 方案二单符号功率=%.2f, 方案一单符号功率=%.2f\n', factor_s2, factor_s1);

% 初始化用于存储结果的cell数组
sinr_results = cell(1, 2);
data_indices_results = cell(1, 2);

%% 2. 循环计算两个方案
for scheme_to_analyze = 1:2
    
    % --- 根据方案设定调制方式、功率和数据排布 ---
    if scheme_to_analyze == 1
        disp('=====================================================');
        disp('开始分析【方案一】，使用 QPSK 调制...');
        modu = 4; % QPSK
        factor = factor_s1; % 使用为方案一计算出的功率
        dd_template = ones(M, N);
        dd_template(1:l_max+1, :) = 0;
        dd_template(M-l_max:M, :) = 0;
    elseif scheme_to_analyze==2
        disp('=====================================================');
        disp('开始分析【方案二】，使用 BPSK 调制...');
        modu = 2; % BPSK
        factor = factor_s2; % 使用基准功率
        dd_template = ones(M, N);
        dd_template(1:l_max+1, 1:2*k_max+1) = 0;
        dd_template(1:l_max+1, N-2*k_max+1:N) = 0;
        dd_template(M-l_max+1:M, 1:2*k_max+1) = 0;
        dd_template(M-l_max+1:M, N-2*k_max+1:N) = 0;          
    end
    
    % --- 生成单次快照的实际发射信号 ---
    bit_stream = randi([0, modu-1], M, N);
    x_modulated = qammod(bit_stream, modu, 'UnitAveragePower', true);
    x_modulated = x_modulated * sqrt(factor); % 使用对应方案的功率定标
    x_modulated(dd_template == 0) = 0;
    dd_vec = x_modulated(:);
    
    data_indices = find(dd_template(:) ~= 0);
    data_indices_results{scheme_to_analyze} = data_indices;
    
    % --- 核心计算 (构建H_eff, G, GH) ---
    disp('正在执行核心计算...');
    hw=zeros(M,N);
    for l=0:M-1, for k=0:N-1, for i=1:P, theta=exp(1j*2*pi*ki(i)*(li(i))/(M*N)); delta_term = (l == li(i)); hw(l+1,k+1)=hw(l+1,k+1)+h_p(i)*h_exp(i)*delta_term*zeta_N(k-ki(i),N).*theta; end, end, end
    H_eff = zeros(M*N, M*N);
    for l=0:M-1, for k=0:N-1, row_idx=k*M+(l+1); for l_prime=0:M-1, for k_prime=0:N-1, col_idx=k_prime*M+(l_prime+1); h_l=mod(l-l_prime,M); h_k=mod(k-k_prime,N); H_eff(row_idx,col_idx)=hw(h_l+1,h_k+1); end, end, end, end
    heff = H_eff'; Ie = eye(size(H_eff,2));

    % !!! 关键修正: MMSE均衡器使用对应方案的factor !!!
    % a = 噪声功率 / 信号功率。原始噪声功率为1，信号功率为factor。
    G = (heff * ((H_eff * heff + (1/factor) * Ie))^(-1));
    
    GH = G * H_eff;
    GG_h = G * G';
    disp('核心矩阵计算完成。');
    
    % --- 计算瞬时S, I, N以及SINR ---
    disp('正在计算瞬时SINR...');
    total_received_vec = GH * dd_vec;

    signal_only_vec = diag(diag(GH)) * dd_vec; 
    interference_vec = total_received_vec - signal_only_vec;
    
    signal_power_inst = abs(signal_only_vec).^2;
    interference_power_inst = abs(interference_vec).^2;
    noise_power_inst = diag(GG_h);
    
    sinr_inst = signal_power_inst ./ (interference_power_inst + noise_power_inst + 1e-20);
    
    sinr_results{scheme_to_analyze} = sinr_inst;
    disp(['方案', num2str(scheme_to_analyze), ' 计算完成。']);
end

%% 3. 结果对比与可视化
% --- 提取两个方案对应数据符号的SINR ---
sinr_s1 = sinr_results{1}(data_indices_results{1});
sinr_s2 = sinr_results{2}(data_indices_results{2});
% nonzero_sinr_s1 = sinr_s1(sinr_s1 ~= 0);
% nonzero_sinr_s2 = sinr_s2(sinr_s2 ~= 0);
nonzero_sinr_s1 = sinr_s1(sinr_s1 > 0); 
% 计算平均值
avg_sinr_s1_linear = mean(nonzero_sinr_s1);

% 对于方案二
nonzero_sinr_s2 = sinr_s2(sinr_s2 > 0);
avg_sinr_s2_linear = mean(nonzero_sinr_s2);


% --- 打印结果 ---
fprintf('\n================== 结果分析 ==================\n');
fprintf('方案一 (QPSK): 数据符号的平均SINR = %.4f \n', avg_sinr_s1_linear );
fprintf('方案二 (BPSK): 数据符号的平均SINR = %.4f \n', avg_sinr_s2_linear);
fprintf('==============================================\n');





%% 4. zeta_N 辅助函数
function output = zeta_N(k, N)
    tolerance = 1e-10;
    if (abs(k) < tolerance || abs(k-N) < tolerance), output = 1;
    else, exp_term = exp(-1i * pi * (N-1)*k/N); numerator = sin(pi*k);
        if abs(numerator)<tolerance, numerator=0; end; denominator = sin(pi*k/N);
        if abs(denominator)<tolerance, output=0; else, output=exp_term*(numerator/(N*denominator)); end;
    end
end