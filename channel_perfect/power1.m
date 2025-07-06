clc;clear;close all;

M=64;
N=32;
fc=32e9;delta_f=120e3;
u_max=120*1000/3600;c=3e8;
%计算抽头
delay=[30,150,310,370,710,1090,1730,2510]*1e-9;
delay_tap0=delay*M*delta_f;
li=round(delay_tap0);
doppler_max=fc*u_max/c;
theta0= -pi + 2*pi*rand(1, length(delay));
doppler=doppler_max*cos(theta0) ;
ki=doppler*N*(1/delta_f);
ki=round(ki,2);
% EVA信道，时延为0的路径是增益最大的，以后依次递减
h_p_db=[1.5,1.4,3.6,0.6,9.1,7,12,16.9];
h_p = 10.^((-1) * h_p_db/10)*10;%假设第一个路径增益为10
l_max=20;k_max=2;
P=length(delay);%路径数
modu=4;%BPSK
% 方案一的速率
bit= (M-2*l_max-1)*N*log2(modu);
% 方案一的数据个数
ans1=(M-2*l_max-1)*N;

%% 符号检测SNR
SNR_persym=8:2:14;
iter=[4e3,4e4,4e5,4e6];
ans2=(M-2*l_max-1)*N+(N-2*k_max-1)*(2*l_max+1);
power_persym = 10.^(SNR_persym/10); % 每个符号的能量
power0=SNR_persym*ans2;

ber = zeros(1, length(SNR_persym));

for i_snr=1:length(SNR_persym)
    power=power0(i_snr);
    ber_iter = zeros(1, iter(i_snr));
   
    num_frames =ceil(iter(i_snr)/bit);
    frame_errors = zeros(1, num_frames);
    frame_bits = zeros(1, num_frames);
    parfor f=1:num_frames
        disp(f)
        h_exp=exp( 1j*2 * pi * rand(1, P));%生成路径相位
        h_exp_round = round(h_exp, 4);%取小数点后四位
        
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声
        noise=reshape(noise,M,N);

        %% 生成bit流
       
        bit_all = randi(modu,M,N)-ones(M,N); %DD域
        x = qammod(bit_all,modu,'UnitAveragePower',true);%调制后的信号的平均功率为 1，总能量/个数，能量为模长平方
        dd=reshape(x,M,N);              
        dd=dd.*sqrt((power/ans1));       

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
        
        %% 方案1：信道可以完美估计
        for i=1:l_max+1
            for j=1:N
                dd(i,j)=0;
            end
        end
        for i=M-l_max+1:1:M
            for j=1:N
                dd(i,j)=0;
            end
        end
        Xp=sqrt(10^(30/10));    dd(1,1)=Xp;
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
     
        for l = 0:M-1
            for k = 0:N-1
                for l_prime = 0:M-1
                    for k_prime = 0:N-1
                        % 计算 h_w 的索引值（带有周期性边界）
                        h_index_l = mod(l - l_prime, M) ;  % 对 l 进行周期性索引
                        h_index_k = mod(k - k_prime, N) ;  % 对 k 进行周期性索引
                        Y(l+1,k+1) = Y(l+1,k+1) + dd(l_prime+1, k_prime+1).* hw(h_index_l+1, h_index_k+1);
                    end
                end
            end
        end
        yabs0=abs(Y);
        %% 噪声
        Y=Y+noise;
               

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
        %验证是否含等效信道，含噪声的接收还是Y
        y=H_eff*dd_vec;
        yabs1=reshape(abs(y),M,[]);

         %% --- 导频干扰消除 (PIC) ---
        
        % 步骤 1: 重建一个只包含导频的发送矩阵
        % 创建一个全零矩阵，只在导频位置(1,1)放置导频值Xp
        dd_pilot_only = zeros(M, N);
        dd_pilot_only(1, 1) = Xp;

        % 步骤 2: 估计导频部分对接收信号的贡献  
        % 使用高效的FFT方法来完成这个卷积
        Y_pilot_est = ifft2(fft2(dd_pilot_only) .* fft2(hw));

        % 步骤 3: 从总接收信号中减去估计的导频影响
        % Y_data_only 就是消除了导频干扰后，留给数据检测器处理的信号
        Y_data_only = Y - Y_pilot_est;
        
        %% MMSE 检测 (DD域)
        H_lmmse = (H_eff' * H_eff + 1 * eye(size(H_eff, 2))) \ H_eff';
        dd_est_matrix=H_lmmse*Y_data_only(:);
        dd_est_matrix= reshape(dd_est_matrix, M, []);
      
        % --- 步骤 5: 提取数据符号并解调 ---
        % 首先，找到原始数据符号的位置 (不包括值为0的保护带和值为Xp的导频)
        data_indices = find(dd ~= 0 & dd ~= Xp);
        % 从恢复的信号矩阵中，提取出对应数据位置的符号
        dd_est_vec = dd_est_matrix(data_indices);
        dd_est_vec=dd_est_vec./sqrt((power/ans1));
        % 提取原始发送的比特流 (用于后续的误码率比较)
        bit_orig = bit_all(data_indices);
        % 将恢复的符号解调成比特流
        bit_est = qamdemod(dd_est_vec, modu, 'UnitAveragePower', true);
        % --- 步骤 6: 计算当前迭代的比特误码率 (BER) ---
        [num_errors_this_frame, ~] = biterr(bit_orig(:), bit_est(:));
        num_bits_this_frame = length(bit_orig(:));
        frame_errors(f) = num_errors_this_frame;
        frame_bits(f) = num_bits_this_frame;
        

   end
   total_bit_errors = sum(frame_errors);
   total_bits_transmitted = sum(frame_bits);
   ber(i_snr) =total_bit_errors / total_bits_transmitted;% 计算当前信噪比下的平均BER
end

figure;
semilogy(SNR_persym, ber, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
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


