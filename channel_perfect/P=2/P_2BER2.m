clc;clear;close all;

M=64;
N=32;
fc=64e9;delta_f=120e3;
u_max=120*1000/3600;c=3e8;
%计算抽头
delay=[30,150]*1e-9;
delay_tap0=delay*M*delta_f;
li=round(delay_tap0);
doppler_max=fc*u_max/c;
theta0= -pi + 2*pi*rand(1, length(delay));
doppler=doppler_max*cos(theta0) ;
ki=doppler*N*(1/delta_f);
k_max0=doppler_max*N*(1/delta_f);

% EVA信道，时延为0的路径是增益最大的，以后依次递减
h_p_db=[1.5,1.4];
relative_power_linear = 10.^((-1) * h_p_db/10);%db转换为功率
total_power = sum(relative_power_linear);
h_p = relative_power_linear / total_power;%每一个增益所占据的百分比
h_exp=exp( 1j*2 * pi * rand(1, length(h_p_db)));

l_max=20;k_max=2;
P=length(delay);%路径数
modu=2;%BPSK
%方案二的数据个数
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
bit= ans2*log2(modu);


%% 符号检测SNR
SNR=10:2:16;
iter=[4e4,4e4,4e5,4e5];

power_persym = 10.^(SNR/10); % 每个符号的能量，考虑到高阶调制
power=power_persym*ans2;%高阶调制后的全部符号能量
factor=power./ans2;%计算真正每个符号的能量
ber = zeros(1, length(SNR));

for i_snr=1:length(SNR)
    ber_iter = zeros(1, iter(i_snr));%一个frame错误bit数
    num_frames =ceil(iter(i_snr)/bit);
    frame_errors = zeros(1, num_frames);
    frame_bits = zeros(1, num_frames);
    parfor f=1:num_frames
        disp(f)
        %% 生成bit流
        bit_all = randi(modu,M,N)-ones(M,N); %DD域
        x = qammod(bit_all,modu,'UnitAveragePower',true);%调制后的信号的平均功率为 1，总能量/个数，能量为模长平方

        dd=reshape(x,M,N);
        % 方案二的功率分配
        dd=dd.*sqrt(factor(i_snr));

        %保护间隔
        for i=1:l_max+1
            for j=1:2*k_max+1
                dd(i,j)=0;
            end
        end
        for i=1:l_max+1
            for j=N:-1:N-2*k_max+1
                dd(i,j)=0;
            end
        end
        for i=M-l_max+1:1:M
            for j=1:2*k_max+1
                dd(i,j)=0;
            end
        end
        for i=M-l_max+1:M
            for j=N:-1:N-2*k_max+1
                dd(i,j)=0;
            end
        end

        hw=zeros(M,N);Y=zeros(M,N);

        %% 分数多普勒域下，DD域等效信道代码,发送数据和接受数据都是dd域
        for l=0:M-1
            for k=0:N-1
                for i=1:P
                    theta=exp(1j*2*pi*ki(i)*(li(i))/(M*N));
                    delta_term = (l == li(i)); % 如果 l 等于 li(i)，则 delta 为 1
                    hw(l+1,k+1)=hw(l+1,k+1)+h_p(i)*h_exp(i)*delta_term*zeta_N(k-ki(i),N).*theta;
                end
            end
        end
        mag_hw = abs(hw);



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
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声
        y=H_eff*dd_vec+noise;
        % y=H_eff*dd_vec;
        yabs=reshape(abs(y),M,[]);



        %% MMSE 检测
        % H_lmmse = (H_eff' * H_eff + 1 * eye(size(H_eff, 2)))^(-1)* H_eff';
        heff=H_eff';
        Ie=eye(size(H_eff, 2));
        H_lmmse=(heff*((H_eff*heff+(1/factor(i_snr))*Ie))^(-1));
        dd_est_matrix=H_lmmse*y;
        dd_est_matrix= reshape(dd_est_matrix, M, []);

        % --- 步骤 5: 提取数据符号并解调 ---
        % 找到原始数据符号的位置 (排除保护带的0和导频Xp)
        data_indices = find(dd ~= 0 );

      
        % 提取对应位置的恢复符号
        dd_est_vec = dd_est_matrix(data_indices);
        dd_vec=dd(data_indices);

        dd_est_vec=dd_est_vec./sqrt(factor(i_snr));
        % 提取原始发送的比特 (用于比较)
        bit_orig = bit_all(data_indices);

        % 解调恢复的符号，得到估计的比特
        bit_est = qamdemod(dd_est_vec, modu, 'UnitAveragePower', true);

        % --- 步骤 6: 计算本次迭代的比特误码率 (BER) ---
        [num_errors_this_frame, ~] = biterr(bit_orig(:), bit_est(:));
        num_bits_this_frame = length(bit_orig(:));
        frame_errors(f) = num_errors_this_frame;
        frame_bits(f) = num_bits_this_frame;
        

        % % 2. 定义“中心数据块”区域并找到其索引
        % data_mask = (dd ~= 0);
        % center_block_mask = false(M, N);
        % center_block_mask(l_max+2 : M-l_max, :) = true;
        % center_indices = find(center_block_mask & data_mask);
        % 
        % % 3. 定义“边缘数据块”区域 (即顶部和底部的数据)并找到其索引
        % % 注意：中心块的掩码取反，再和总数据掩码相与，即可得到所有边缘数据
        % edge_indices = find(~center_block_mask & data_mask);
        % 
        % % 5. 提取并解调“中心块”数据
        % dd_est_center = dd_est_matrix(center_indices);
        % bit_orig_center = bit_all(center_indices);
        % bit_est_center = qamdemod(dd_est_center./sqrt(factor(i_snr)), modu, 'UnitAveragePower', true);
        % [num_errors_center, num_bits_center] = biterr(bit_orig_center(:), bit_est_center(:));
        % frame_errors_center(f) = num_errors_center;
        % frame_bits_center(f) = num_bits_center;   
        % num_bits_center
        % 
        % % 6. 提取并解调“边缘块”数据
        % dd_est_edge = dd_est_matrix(edge_indices);
        % bit_orig_edge = bit_all(edge_indices);
        % bit_est_edge = qamdemod(dd_est_edge./sqrt(factor(i_snr)), modu, 'UnitAveragePower', true);
        % [num_errors_edge, num_bits_edge] = biterr(bit_orig_edge(:), bit_est_edge(:));
        % frame_errors_edge(f) = num_errors_edge;
        % frame_bits_edge(f) = num_bits_edge;
        % num_bits_edge
        
    end
    total_bit_errors = sum(frame_errors);
    total_bits_transmitted = sum(frame_bits);
    ber(i_snr) =total_bit_errors / total_bits_transmitted;% 计算当前信噪比下的平均BER
end

figure;
semilogy(SNR, ber, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('信噪比 (dB)');ylabel('比特误码率 (BER)');title('方案二MMSE检测后系统BER性能曲线');legend('MMSE-BER');

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


