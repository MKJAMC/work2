clc;clear;close all;
%% 方案二
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
modu=2;%BPSK
%方案二的数据个数
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
bit= ans2*log2(modu);


%% 符号检测SNR
SNR=10:2:18;
iter=[4e4,4e4,4e5,4e5,4e6];
% SNR=18;iter=[4e5];
chanel_err=zeros(1,length(SNR));
power_persym = 10.^(SNR/10); % 每个符号的能量，考虑到高阶调制
power=power_persym*ans2;%高阶调制后的全部符号能量
factor=power./ans2;%计算真正每个符号的能量
ber = zeros(1, length(SNR));
xp=sqrt(power_persym(1)*1e5);%30dB是1e3

for i_snr=1:length(SNR)
    ber_iter = zeros(1, iter(i_snr));%一个frame错误bit数
    num_frames =ceil(iter(i_snr)/bit);
    frame_errors = zeros(1, num_frames);
    frame_bits = zeros(1, num_frames);
    parfor f=1:num_frames
        disp(f)
        theta0= -pi + 2*pi*rand(1, length(delay));
        doppler=doppler_max*cos(theta0) ;
        ki=doppler*N*(1/delta_f);
        h_exp=exp( 1j*2 * pi * rand(1, length(h_p_db)));

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
        dd(1,1)=xp;
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

        % 频域相乘
        H_freq = fft2(hw);
        % 步骤 B: 对输入 dd 做二维傅里叶变换
        dd_freq = fft2(dd);
        % 步骤 C: 在频域进行逐元素相乘 (这步等效于 H_eff * dd_vec)
        y_freq = H_freq .* dd_freq;
        % 步骤 D: 将结果通过二维逆傅里叶变换返回到原始域
        y = ifft2(y_freq);
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声
        noise=reshape(noise,M,N);
        Y=y+noise;
        YABS=abs(Y);
        [max_values0, max_indices0] = max(YABS, [], 2);%选出了每一行的最大值

         %% 互相关信道估计
        step=0.1;
        yd=zeros(1,size(Y,2));%待抽取出的每一行
        colnum=length(-k_max:step:k_max);%step细化以后的总列数
        r=zeros(l_max+1,colnum);%相关值进行存放，行代表着带估计的时延位置，列代表着带估计的分数多普勒的位置
        %计算每一条路径的相关位置
        for l=1:l_max+1
            jj=1; % <-- 我们定义了一个独立的计数器变量 jj
            yd=Y(l,:);
            for stepval=-k_max:step:k_max  % <-- 'stepval' 是这个循环的循环变量
                for k=1:size(Y,2)
                    r(l,jj)=r(l,jj)+yd(k)*conj(zeta_N(k-stepval-1,N));
                end
                jj=jj+1; % <-- 我们修改的是独立的计数器 jj，而不是循环变量 stepval
            end
        end
        rabs=abs(r);
        [max_values, max_indices] = max(rabs, [], 2);%选出了每一行的最大值
        li_est = []; h_est=[];h_phi_est=[];ki_est=[];k_cor=[];temp=[];
        [sorted_values, sorted_indices] = sort(max_values0(1:l_max+1), 'descend');
        top_indices = sorted_indices(1:8);
        li_est = sort(top_indices - 1);

        for i=1:length(li_est)
            yd= Y(li_est(i)+1,:); % 选取正确的接收信号行
            k_hat=(-k_max)+step*(max_indices(li_est(i)+1)-1);
            k_cor=[k_cor,k_hat];
            peak_idx = max_indices(li_est(i)+1); % 粗略搜索找到的峰值索引
            % 确定边界索引，并处理边界 保证索引不越界 (例如峰值在第一个或最后一个点)
            left_idx = max(1, peak_idx - 1);
            right_idx = min(colnum, peak_idx + 1);
            % 使用峰值点和它两边的点作为 GSS 的初始边界
            b_l = (-k_max) + step * (left_idx - 1);
            b_u = (-k_max) + step * (right_idx - 1);
            max_iter=50;
            eta = (sqrt(5)-1)/2;
            % 在循环外计算初始的 b1 和 b2 以及它们的函数值
            b1 = b_u - eta * (b_u - b_l);
            b2 = b_l + eta * (b_u - b_l);
            f1_val = 0;
            for k=1:size(Y,2)
                f1_val = f1_val + yd(k)*conj(zeta_N(k-b1-1,N));
            end
            f1_val = abs(f1_val); % 取模

            f2_val = 0;
            for k=1:size(Y,2)
                f2_val = f2_val + yd(k)*conj(zeta_N(k-b2-1,N));
            end
            f2_val = abs(f2_val); % 取模

            for iter_gss = 1:max_iter
                if f1_val > f2_val
                    b_u = b2;
                    b2 = b1;
                    f2_val = f1_val; % 无需重新计算 f2
                    b1 = b_u - eta * (b_u - b_l);
                    f1_val = 0;
                    for k=1:size(Y,2)
                        f1_val = f1_val + yd(k)*conj(zeta_N(k-b1-1,N));
                    end
                    f1_val = abs(f1_val); % 只需要计算一个新的 f1
                else
                    b_l = b1;
                    b1 = b2;
                    f1_val = f2_val; % 无需重新计算 f1
                    b2 = b_l + eta * (b_u - b_l);
                    f2_val = 0;
                    for k=1:size(Y,2)
                        f2_val = f2_val + yd(k)*conj(zeta_N(k-b2-1,N));
                    end
                    f2_val = abs(f2_val); % 只需要计算一个新的 f2
                end
            end
            nu_opt = (b_l + b_u) / 2;
            ki_est=[ki_est,nu_opt];
            yd = Y(li_est(i)+1, :); % 获取当前路径所在行的接收信号
            r_refined = 0;
            for k=1:N
                r_refined = r_refined + yd(k) * conj(zeta_N(k-nu_opt-1, N));
            end
            h_est_i = abs(r_refined) / (xp); % 注意这里是 xp 而不是 power_persym
            h_est = [h_est, h_est_i];

            % 2. 估计相位 (使用 r_refined 的相位)
            if abs(r_refined) > 1e-9 % 避免除以0
                h_phi_est_i = r_refined / abs(r_refined);
            else
                h_phi_est_i = 1;
            end
            h_phi_est = [h_phi_est, h_phi_est_i];
        end

        %% 估计的数值生成信道,NMSE
        hw_est=zeros(M,N);
        for l=0:M-1
            for k=0:N-1
                for i=1:P
                    theta=exp(1j*2*pi*ki_est(i)*(li_est(i))/(M*N));
                    delta_term = (l == li_est(i)); % 如果 l 等于 li(i)，则 delta 为 1,l是从0开始计算
                    hw_est(l+1,k+1)=hw_est(l+1,k+1)+h_est(i)*h_phi_est(i)*delta_term*zeta_N(k-ki_est(i),N).*theta;
                end
            end
        end
        %% MMSE 检测
       
        y_freq = fft2(Y); % y已经是矩阵形式，无需reshape
        alpha = 1 / factor(i_snr);
        H_freq = fft2(hw_est);
        H_mmse_freq = conj(H_freq) ./ (abs(H_freq).^2 + alpha);
        z_freq_initial = H_mmse_freq .* y_freq;
        dd_est_matrix = ifft2(z_freq_initial);

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


