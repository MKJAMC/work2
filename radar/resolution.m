clc; clear; close all;

%% ================= 第一部分：参数设置 =================
M = 64;
N = 32;
fc = 64e9;          
delta_f = 120e3;    
c = 3e8;            

% 基础参数
delay_res = 1 / (M * delta_f);      
doppler_res = 1 / (N * (1/delta_f));          
distance_res = c * delay_res / 2;   
v_res = doppler_res * c / 2 / fc;   
R_max = M * distance_res;           
V_max = N * v_res / 2;   %存在正负速度           

% 样本生成参数
datanums = 4600;    % 总样本数
num_targets = 3;    % 目标数量

% 初始化存储
Dataset_Y = zeros(M, N, 2, datanums); % [M, N, 2, Batch]
Target_Delay_Index = zeros(datanums, num_targets);
Target_Doppler_Index = zeros(datanums, num_targets);

fprintf('1. 正在生成目标参数...\n');
% 模拟生成目标参数 (分数索引)
for i = 1:datanums
    % 随机生成 0.1 ~ 0.9 倍范围内的目标
    Target_Range = sort(rand(1, num_targets) * R_max * 0.9);
    Target_Velocity = sort((rand(1, num_targets)*2-1) * V_max * 0.8);
    
    % 转换为分数索引 (Fractional Indices)
    % li: 对应文档中的 \iota_{\tau_i} (分数)
    % ki: 对应文档中的 \kappa_{\nu_i} (分数)
    Target_Delay_Index(i, :) = Target_Range ./ distance_res;
    Target_Doppler_Index(i, :) = Target_Velocity ./ v_res;
    
end
%% ================= 第二部分：OTFS 信号生成 (基于文档公式) =================
fprintf('2. 正在生成 OTFS 回波信号 Y (基于文档公式 #10-#13)...\n');

% --- OTFS 调制参数 ---
l_max = 20; k_max = 2;
cen_modu = 2; guard_modu = 4;
bits_per_qam_symbol = log2(cen_modu);
bits_per_guard_symbol = log2(guard_modu);
g = 4; 
SNR_range = 10:2:18; 

% 预计算部分参数
mo = mod((N-4*k_max-1), g);
q1 = floor((N-4*k_max-1)/ g);
guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
data_rows = (l_max + 2):(M - l_max);
guard_rows = [1:(l_max + 1), (M - l_max + 1):M];

% 开始循环生成
% parfor i_data = 1:datanums % 如需加速可开启并行
for i_data = 1:datanums
    if mod(i_data, 100) == 0
        fprintf('处理进度: %d / %d\n', i_data, datanums);
    end
    
    % --- A. 获取当前样本的目标参数 ---
    li = Target_Delay_Index(i_data, :); % 1x3, 分数时延索引 \iota
    ki = Target_Doppler_Index(i_data, :); % 1x3, 分数多普勒索引 \kappa
    
    % --- B. 生成信道增益 (瑞利分布) ---
    % 对应公式 #6 中的 h_i e^{j\phi_i}
    % 实部虚部服从高斯分布 -> 模服从瑞利分布
    h_complex = (randn(1, num_targets) + 1j * randn(1, num_targets)) / sqrt(2);
    
    % --- C. 生成发射信号矩阵 dd (X) ---
    dd = zeros(M, N);
    
    % C1. 中央数据区 (IM + QAM)
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
    
    % C2. 保护间隔区
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
    
    % C3. 功率与导频
    current_SNR = SNR_range(randi(length(SNR_range)));
    Target_SNR(i_data) = current_SNR;
    xp = 1000; %xp为固定数值
    dd(1,1) = xp; % 放置导频
    
    % --- D. 生成 DD 域等效信道 hw ---
    % 依据文档公式 #10, #11, #12, #13
    % hw[l,k] 表示信道冲激响应在 DD 网格上的分布
    hw = zeros(M, N);
    
    % 网格坐标矩阵 (用于向量化计算 zeta 函数)
    [L_grid, K_grid] = ndgrid(0:M-1, 0:N-1);
    
    for p = 1:num_targets
        tau_idx = li(p);  % \iota_i
        nu_idx = ki(p);   % \kappa_i
        gain = h_complex(p);
        
        % 1. 计算多普勒域扩散 G (公式 #11)
        % 输入为 (k - \kappa_i)，沿着 N 轴
        G_term = zeta_N(K_grid - nu_idx, N);
        
        % 2. 计算时延域扩散 F (公式 #12)
        % 输入为 (l - \iota_i)，沿着 M 轴 (此处体现分数时延的 Sinc 形状)
        F_term = zeta_N(L_grid - tau_idx, M);
        
        % 3. 计算相位偏移 \theta (公式 #13)
        % theta = exp(-j * 2pi * \iota * \kappa / MN)
        theta_term = exp(-1j * 2 * pi * tau_idx * nu_idx / (M*N));
        
        % 4. 累加信道响应
        % h_w = sum( h_i * e^{phi} * G * F * theta )
        hw = hw + gain * G_term .* F_term * theta_term;
    end
    
    % --- E. 通过信道 (频域乘法实现 2D 循环卷积) ---
    % Y_DD = X_DD * H_eff
    H_freq = fft2(hw);
    dd_freq = fft2(dd);
    Y_clean = ifft2(H_freq .* dd_freq);
    
    % --- F. 添加噪声 ---
    sig_power = mean(abs(Y_clean(:)).^2);
    noise_power = sig_power / (10^(current_SNR/10));
    noise = sqrt(noise_power/2) * (randn(M, N) + 1j * randn(M, N));
    
    Y = Y_clean + noise;
    
    % --- G. 存入数据集 ---
    Dataset_Y(:, :, 1, i_data) = real(Y);
    Dataset_Y(:, :, 2, i_data) = imag(Y);

    
%% 绘制DD域分数时延和多普勒的图像
% for p = 1:num_targets
%     fprintf('  目标 %d: 时延索引(li) = %.4f, 多普勒索引(ki) = %.4f\n', ...
%         p, li(p), ki(p));
% end
% figure(1);
% H1 = abs(hw);
% h = bar3(H1); % 创建三维条形图
% 
% colormap(parula); 
% 
% for k = 1:length(h)
%     % 1. 按高度着色（使颜色与幅度对应，而不是按列）
%     zdata = h(k).ZData;
%     h(k).CData = zdata;
%     h(k).FaceColor = 'interp'; % 平滑着色，消除块状感
% 
%     % 2. 解决“灰暗”的关键：弱化边框
%     % 将边框设为非常细的深灰色，或者干脆设为 'none'
%     h(k).EdgeColor = [0.4 0.4 0.4]; 
%     h(k).LineWidth = 0.1; 
% end
% 
% % --- 提升质感：添加环境光 ---
% grid on;
% set(gca, 'GridAlpha', 0.2); % 让网格线淡一点，不干扰视线
% view(-37.5, 30); 
% axis tight;
% 
% % 添加光照，增加表面的金属光泽，去除死板的灰色感
% camlight headlight; 
% lighting gouraud; 
% 
% % 设置背景为纯白，看起来更通透
% set(gcf, 'Color', 'w');
% axis tight;
% xlabel('多普勒 (N)'); ylabel('时延 (M)'); zlabel('幅度');
% title("分数时延和分数多普勒在DD域图像")
end


fprintf('数据生成完成！Dataset_Y 大小: %s\n', mat2str(size(Dataset_Y)));

%% ================= 附录：核心函数 zeta_N =================
function output = zeta_N(k, N)
    % 实现文档中的基函数 G 和 F (Dirichlet Kernel)
    % 公式 #11: 1/N * sum(exp(-j*2pi*n*k/N))
    % 等价于: 1/N * exp(-j*pi*(N-1)*k/N) * sin(pi*k) / sin(pi*k/N)
    
    tolerance = 1e-10;
    
    % 初始化
    output = zeros(size(k));
    
    % 1. 处理 k 为整数倍 N 的峰值情况 (L'Hospital法则极限为1)
    % 利用 mod 判断 k 是否接近 N 的整数倍
    mask_peak = abs(sin(pi * k / N)) < tolerance;
    % 注意：如果 k 是 N 的整数倍，sin(pi*k) 也是 0，极限是 N/N * phase
    % 但这里归一化系数是 1/N，所以幅度是 1
    % 相位项: exp(-j * pi * (N-1) * (integer*N) / N) -> exp(-j * pi * (N-1) * integer)
    % 简单起见，通常 peak 处设为 1 或对应的相位，这里按公式计算：
    
    % 针对峰值的特殊处理 (直接代入公式会有 0/0)
    k_peak = k(mask_peak);
    % 当 k -> 0, output -> 1. 当 k -> N, output -> 1 (周期性)
    % 实际上要保留相位项
    % exp_term_peak = exp(-1i * pi * (N - 1) * k_peak / N);
    % output(mask_peak) = exp_term_peak .* 1; % limit of sin/Nsin is 1
    % 但为了数值稳定，直接设幅值为1的近似：
    output(mask_peak) = 1; 
    
    % 2. 处理非峰值情况 (正常 Sinc 形状)
    mask_normal = ~mask_peak;
    if any(mask_normal(:))
        k_val = k(mask_normal);
        
        % 相位项 (对应文档公式 #11 中的 e^{-j(N-1)\pi \kappa / N})
        exp_term = exp(-1i * pi * (N - 1) * k_val / N);
        
        % 幅度项 (Dirichlet 核)
        numerator = sin(pi * k_val);
        denominator = sin(pi * k_val / N);
        
        output(mask_normal) = exp_term .* (numerator ./ (N * denominator));
    end
end