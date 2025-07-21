clc;clear;close all;
clc;
clear;
close all;

%% 1. 参数定义
M = 64; % DD矩阵的总行数

% 信道延时抽头 (与您代码中计算结果一致)
li = [0, 1, 2, 3, 5, 8, 13, 19];
num_paths = length(li);

% 定义数据区
data_start_row = 22;
data_end_row = 44;

%% 2. 构建发射信号
% 创建一个 M x 1 的向量代表DD矩阵的一列
% 仅在数据区填充1，其余为0 (代表保护间隔)
tx_signal = zeros(M, 1);
tx_signal(data_start_row:data_end_row) = 1;


%% 3. 模拟延时叠加过程
% 初始化一个计数器向量，用于记录每一行上的叠加次数
superposition_counts = zeros(M, 1);

% 遍历所有发射信号为1的行 (即数据区的每一行)
for l_tx = data_start_row:data_end_row
    
    % 让这一行的信号经过所有的延时路径
    for i = 1:num_paths
        
        current_delay = li(i);
        
        % 计算延时后的接收行号 (l_rx)
        % 注意：MATLAB是1-based索引, mod运算是0-based，需转换
        % a. 将发射行号转为0-based (1->0, 2->1, ...)
        l_tx_zero_based = l_tx - 1;
        
        % b. 计算循环移位后的0-based接收行号
        l_rx_zero_based = mod(l_tx_zero_based + current_delay, M);
        
        % c. 转回1-based的MATLAB索引
        l_rx = l_rx_zero_based + 1;
        
        % d. 在对应的接收行上，计数加一
        superposition_counts(l_rx) = superposition_counts(l_rx) + 1;
        
    end
end

%% 4. 显示结果
% 创建一个表格来清晰地展示结果
row_numbers = (1:M)';
results_table = table(row_numbers, superposition_counts, 'VariableNames', {'接收行号', '叠加路径总数'});

% 在命令窗口中显示完整的表格
disp('每一接收行上的有效路径叠加总数:');
disp(results_table);


%% 5. 可视化结果
figure;
bar(row_numbers, superposition_counts, 'FaceColor', '#0072BD');
grid on;
title('每一接收行上的有效路径叠加数量可视化');
xlabel('接收行号 (l_rx)');
ylabel('有效叠加路径数');
xlim([0, M + 1]); % 设置X轴范围
ylim([0, max(superposition_counts) + 1]); % 设置Y轴范围
set(gca, 'FontSize', 12); % 设置坐标轴字体大小

% M=16;
% N=16;
% hw=zeros(M,N);
% P=1;
% dd=zeros(M,N);
% dd(1,:)=10;
% 
% Y=zeros(M,N);
% li=0;ki=0.1;
% h_p_round=1;
% 
% % g=4;%IM调制，每行分成几组
% % dd(1,:)=10;
% % for j = 1:N/g
% %         % 计算当前组的起始列和结束列
% %         start_col = (j-1) * (N/g) + 1;
% %         end_col = j * (N/g);
% %         % 每组只有一个激活
% %         zero_indices = start_col:end_col;
% %         num_zeros = 3;
% %         zero_positions = randperm(length(zero_indices), num_zeros);
% %         % 将随机选择的部分置为 0
% %         dd(1, zero_indices(zero_positions)) = 0;
% % end
% 
% 
% for l=0:M-1
%     for k=0:N-1
%         for i=1:P
%             theta=exp(1j*2*pi*ki(i)*(li(i))/(M*N));
%             delta_term = (l == li(i)); % 如果 l 等于 li(i)，则 delta 为 1
%             hw(l+1,k+1)=hw(l+1,k+1)+h_p_round(i)*delta_term*zeta_N(k-ki(i),N).*theta;
%         end
%     end
% end
% 
%        y_vec = reshape(Y, M*N, 1);
%        dd_vec = reshape(dd, M*N, []);
%         H_eff = zeros(M*N, M*N);
%         for l = 0:M-1       % 遍历矩阵的行
%             for k = 0:N-1   % 遍历矩阵的列
%                 % 计算 Y(l+1, k+1) 在列主序向量 y_vec 中的索引
%                 row_idx = k*M + (l+1);
%                 for l_prime = 0:M-1
%                     for k_prime = 0:N-1
%                         % 计算 dd(l_prime+1, k_prime+1) 在列主序向量 dd_vec 中的索引
%                         col_idx = k_prime*M + (l_prime+1);
%                         h_l = mod(l - l_prime, M);
%                         h_k = mod(k - k_prime, N);
%                         H_eff(row_idx, col_idx) = hw(h_l + 1, h_k + 1);
%                     end
%                 end
%             end
%         end
%         noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声
%         noise=reshape(noise,M,N);
%         y=H_eff*dd_vec;
%         Y=reshape(y,M,[])+noise;
%         H_lmmse = (H_eff' * H_eff + 1 * eye(size(H_eff, 2))) \ H_eff';
%         dd_est_matrix=H_lmmse*Y(:);
%         dd_est_matrix= reshape(dd_est_matrix, M, []);
%         b=abs(dd_est_matrix);
% 
%         figure();
%         plot(1:16,abs(dd(1,:)),'-o');
%         hold on;
%         grid on;
%         plot(1:16,b(1,:),'-^');
%         legend('发射信号','均衡的接收信号');ylabel('幅值');
%         xlabel('满铺');
% 
% function output = zeta_N(k, N)
% %这个函数属于输入一个k，得到一个数值
% tolerance = 1e-10;  % 可调容差
% if (k == N||k==0)
%     output = 1;
% else
%     exp_term = exp(-1i * pi * (N - 1) * k / N);
%     numerator = sin(pi * k);  % 分子 sin(pi * k)
%     if abs(numerator) < tolerance
%         numerator = 0;
%     end
%     denominator = sin(pi * k / N);  % 分母 sin(pi * k / N)
%     output = exp_term * (numerator / (N * denominator));
% end
% end
% 


% Parameters
% N = 8; % Example value, modify as needed
% M = 8; % Example value, modify as needed
% T = 0.1; % Pulse duration, example value
% Delta_f = 1/T; % Frequency difference, example value
%
%
% % t = 0:T:(N-1)*T; % Example time range, modify as needed
% t=0:0.01:0.8;
% % Initialize s(t) signal
% s_t = zeros(size(t));
% s_t_hun=zeros(size(t));
% bit_all = randi(4,M,N)-ones(M,N); %DD域
% modu=4;
% x = qammod(bit_all,modu,'UnitAveragePower',true);%星座图包括的点
% x=zeros(M,N);
% x(2,1)=1;
% m = size(x,1);
% n = size(x,2);
% X_TF =ifft((fft(x)).').* sqrt(m/n); %SFFT -> IDFT in delay domain and DFT in Doppler domain
% % X_T=reshape(X_T,[],1);
% Y = ifft(x, N, 2);


