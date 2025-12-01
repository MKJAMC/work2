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
modu=4;%QPSK
% 方案一的数据个数
ans1=(M-2*l_max-1)*N;
% 方案一的速率
bit= (M-2*l_max-1)*N*log2(modu);
%方案二的数据个数,其MN,l_max,k_max SendSocke不变则个数不变
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
%% 信道估计SNR
SNR = 10:2:18;
iter=[2000,2000,2000,2000,2e3];%迭代次数

chanel_err=zeros(1,length(SNR));
% 方案一的数据个数
ans1=(M-2*l_max-1)*N;
% 方案一的速率
bit= (M-2*l_max-1)*N*log2(modu);
%方案二的数据个数,其MN,l_max,k_max不变则个数不变
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
power_persym=10.^(SNR/10);%每个符号的功率
power=power_persym*ans2;%ans2为基准，求得总功率，方案二的功率为1
powerdb=10*log10(power)

factor=power./ans1;%计算真正每个符号的能量
xp=sqrt(power_persym(1)*1e5);

for i_snr=1:length(SNR)
    nmse=zeros(1,iter(i_snr));
    nums=0;
    for nums=1:iter(i_snr)%每次都重新生成信道
        theta0= -pi + 2*pi*rand(1, length(delay));
        doppler=doppler_max*cos(theta0) ;
        ki=doppler*N*(1/delta_f);

        % ki=[  0.56477 ,    -0.14871  ,   0.28003,0.41167 ,    0.49664 ,     1.7472 ,   -0.82135 ,   -0.77633 ];
        h_exp=exp( 1j*2 * pi * rand(1, length(h_p_db)));
        % h_p = ones(1, P) / P;
        disp(nums)
        %% 生成bit流
        bit_all = randi(modu,M,N)-ones(M,N); %DD域
        x = qammod(bit_all,modu,'UnitAveragePower',true);%调制后的信号的平均功率为 1，总能量/个数，能量为模长平方
        dd=reshape(x,M,N);
        dd=dd.*sqrt((factor(i_snr)));
        %保护间隔
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
        dd(1,1)=xp;
        hw=zeros(M,N);Y=zeros(M,N);

        %% 分数多普勒域下，DD域等效信道代码,发送数据和接受数据都是dd域
        for l=0:M-1
            for k=0:N-1
                for i=1:P
                    theta=exp(1j*2*pi*ki(i)*(li(i))/(M*N));
                    delta_term = (l == li(i)); % 如果 l 等于 li(i)，则 delta 为 1
                    hw(l+1,k+1)=hw(l+1,k+1)+h_p(i)*h_exp(i)*delta_term*zeta_N(k-ki(i),N).*theta;
                    % hw(l+1,k+1)=hw(l+1,k+1)+delta_term*zeta_N(k-ki(i),N).*theta;
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
        %% 噪声
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声
        noise=reshape(noise,M,N);        noiseabs=abs(noise);
        Y=y+noise;        YABS=abs(Y);
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

        %% 黄金分割法
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
            % nu_opt=round(nu_opt,1);
            ki_est=[ki_est,nu_opt];
            % ze=zeros(1,N);
            % for j=1:N
            %     ze(j)=zeta_N(j-1-nu_opt,N);%已知分数多普勒的时候多普勒维度的扩展
            % end
            % [maxze1,ze1index]=max(ze);
            % h_est_i=abs(max(Y(i,:)))/abs((maxze1))/xp(i_snr);
            % h_est=[h_est,h_est_i];
            % temp=r(i, max_indices(i));
            % h_phi_est_i=temp/abs(temp);
            % h_phi_est=[h_phi_est,h_phi_est_i];
            yd = Y(li_est(i)+1, :); % 获取当前路径所在行的接收信号
            r_refined = 0;
            for k=1:N
                r_refined = r_refined + yd(k) * conj(zeta_N(k-nu_opt-1, N));
            end

            % --- 步骤B: 基于这个“干净”的复数相关值来估计增益和相位 ---

            % 1. 估计增益/幅度 (使用 r_refined 的幅度)
            % 理论上，相关峰值的幅度约等于 |h| * xp
            % N 是相关操作的长度，也是相干累积带来的增汁
            h_est_i = abs(r_refined) / (xp); % 注意这里是 xp 而不是 power_persym
            h_est = [h_est, h_est_i];

            % 2. 估计相位 (使用 r_refined 的相位)
            if abs(r_refined) > 1e-9 % 避免除以0
                theta_est = exp(1j * 2 * pi * nu_opt * li_est(i) / (M*N));
                h_phi_est_i = r_refined / abs(r_refined)/theta_est;
            else
                h_phi_est_i = 1;
            end
            h_phi_est = [h_phi_est, h_phi_est_i];
        end

        %
        %
        % % for i=1:length(li_est)
        % %     % % --- 抛物线插值法部分---
        % %     peak_idx = max_indices(li_est(i)+1);
        % %     k_hat=(-k_max) + step*(peak_idx-1);
        % %
        % %     if peak_idx > 1 && peak_idx < colnum
        % %         y_center = rabs(li(i)+1, peak_idx);
        % %         y_left   = rabs(li(i)+1, peak_idx - 1);
        % %         y_right  = rabs(li(i)+1, peak_idx + 1);
        % %         denominator = 2*y_center - y_left - y_right;
        % %         if abs(denominator) > 1e-6
        % %             p = (y_right - y_left) / (2 * denominator);
        % %             nu_opt = k_hat + p * step;
        % %         else
        % %             nu_opt = k_hat;
        % %         end
        % %     else
        % %         nu_opt = k_hat;
        % %     end
        % %     ki_est=[ki_est,nu_opt];
        % %     k_cor=[k_cor,k_hat]; % 保留您的粗略估计用于对比
        % %
        % %      ze=zeros(1,N);
        % %      for j=1:N
        % %          ze(j)=zeta_N(j-1-k_hat,N);%已知分数多普勒的时候多普勒维度的扩展
        % %      end
        % %      [maxze1,ze1index]=max(ze);
        % %      h_est_i=abs(max(Y(i,:)))/abs((maxze1))/xp(i_snr);
        % %      h_est=[h_est,h_est_i];
        % %      temp=r(i, max_indices(i));
        % %      h_phi_est_i=temp/abs(temp);
        % %      h_phi_est=[h_phi_est,h_phi_est_i];
        % %
        % % end
        % % abs(ki_est-ki)
        % % abs(k_cor-ki)
        disp(['不同路径的多普勒=',num2str(ki)]);
        disp(['估计的多普勒抽头=',num2str(ki_est)]);
        disp(['估计的小数点后一位多普勒抽头=',num2str(k_cor)]);
        % disp(['不同路径的时延=',num2str(li)]);
        % disp(['估计的时延抽头=',num2str(li_est)]);
        disp(['路径的信道增益=',num2str(h_p)]);
        disp(['估计的信道增益=',num2str(h_est)]);
        disp(['路径的信道相位=',num2str(h_exp)]);
        disp(['估计的信道相位=',num2str(h_phi_est)]);


        %% 估计的数值生成信道,NMSE
        hw_est=zeros(M,N);
        for l=0:M-1
            for k=0:N-1
                for i=1:length(li_est)
                    theta=exp(1j*2*pi*ki_est(i)*(li_est(i))/(M*N));
                    delta_term = (l == li_est(i)); % 如果 l 等于 li(i)，则 delta 为 1,l是从0开始计算
                    hw_est(l+1,k+1)=hw_est(l+1,k+1)+h_est(i)*h_phi_est(i)*delta_term*zeta_N(k-ki_est(i),N).*theta;
                    % hw_est(l+1,k+1)=hw_est(l+1,k+1)+delta_term*zeta_N(k-ki_est(i),N).*theta;
                end
            end
        end
        num = norm(hw_est - hw, 'fro')^2;    % Σ_{l,k} |hw_est(l,k)-hw(l,k)|^2
        den = norm(hw, 'fro')^2;         % Σ_{l,k} |hw(l,k)|^2
        nmse(nums) = num / den;
    end
    chanel_err(i_snr)=mean(nmse);
end


figure;
semilogy(SNR, chanel_err, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('信噪比 (dB)');ylabel('NMSE');title('NMSE性能曲线');legend('方案1');
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



