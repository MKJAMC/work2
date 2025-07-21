clc;clear;close all;

M=16;
N=16;
l_max=2;k_max=2;
P=2;%路径数
Xp=10;
modu=2;%BPSK
%方案二的数据个数
ans2=(M-2*l_max-1)*N+(N-2*k_max-1)*(2*l_max+1);
bit= ans2*log2(modu);
ans1=(M-2*l_max-1)*N;

%% 信道检测
SNR_dB = 0:5:25;
iter=[2000,2000,2000,2000,2000,2000];%迭代次数

%% 符号检测SNR
% SNR_dB=20;
% iter=[4e6];

ber = zeros(1, length(SNR_dB));
chanel_err=zeros(1,length(SNR_dB));

for i_snr=1:length(SNR_dB)
    SNR = 10^(SNR_dB(i_snr)/10);
    nmse=zeros(1,iter(i_snr));ber_iter = zeros(1, iter(i_snr));
    nums=0;
    total_bits_transmitted = 0;%传输bit
    total_bit_errors = 0;%错误bit
    for nums=1:iter(i_snr)%每次都重新生成信道
        nums
    % while  (nums*bit )<iter(i_snr)
    %     nums = nums + 1; % 迭代次数加1
        li = randperm(l_max, P)';
        k = randi([-k_max,k_max],P,1);kv=randi([-5,5],P,1)*0.1;ki=k+kv;
        % EVA信道，时延为0的路径是增益最大的，以后依次递减
        h_p_db=[1.5,1.4];
        h_p = 10.^((-1) * h_p_db/10)*10;%假设第一个路径增益为10
        h_p_round     = round(h_p, 4);
        h_exp=exp( 1j*2 * pi * rand(1, P));%生成路径相位
        h_exp_round = round(h_exp, 4);
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声
        noise=reshape(noise,M,N);


        %% 生成bit流

        bit_all = randi(modu,M,N)-ones(M,N); %DD域
        x = qammod(bit_all,modu,'UnitAveragePower',true);%调制后的信号的平均功率为 1，总能量/个数，能量为模长平方
        dd=reshape(x,M,N);


        % 方案二的功率消减
        dd=dd.*sqrt(SNR/ans1);
        
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


        Xp=10;    dd(1,1)=Xp;
        hw=zeros(M,N);Y=zeros(M,N);

        %% 分数多普勒域下，DD域等效信道代码,发送数据和接受数据都是dd域
        for l=0:M-1
            for k=0:N-1
                for i=1:P
                    theta=exp(1j*2*pi*ki(i)*(li(i))/(M*N));
                    delta_term = (l == li(i)); % 如果 l 等于 li(i)，则 delta 为 1
                    hw(l+1,k+1)=hw(l+1,k+1)+h_p_round(i)*h_exp_round(i)*delta_term*zeta_N(k-ki(i),N).*theta;
                end
            end
        end
        mag_hw = abs(hw);
        % figure();
        % plot(0:M-1,magnitude(3,:),'bo');hold on;grid on;

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
        yabs=abs(Y);
        Y=Y+noise;

        step=0.1;
        yd=zeros(1,size(Y,2));%待抽取出的每一行
        colnum=((2*k_max+1)/(step))+1;%step细化以后的总列数
        r=zeros(l_max+1,colnum);%相关值进行存放，行代表着带估计的时延位置，列代表着带估计的分数多普勒的位置
        for l=1:l_max+1
            j=1;
            yd=Y(l,:);
            for stepval=-k_max-0.5:step:k_max+0.5
                for k=1:size(Y,2)
                    r(l,j)=r(l,j)+yd(k)*conj(zeta_N(k-stepval-1,N));   %计算出不同step下的自相关值
                end
                j=j+1;
            end
        end

        li_est = []; h_est=[];h_phi_est=[];ki_est=[];
        rabs=abs(r);
        [max_values, max_indices] = max(rabs, [], 2);%选出了每一行的最大值
        for i=2:length(max_values)%假设信道时延不为0
            if(max_values(i)>3)%判断时延是否存在的证据
                li_est=[li_est,i-1];
                ki_est_i=(-k_max-0.5)+step*(max_indices(i)-1);%确定时延以后，得到对应时延的多普勒数值
                ki_est=[ki_est,ki_est_i];
                ze=zeros(1,N);
                for j=1:N
                    ze(j)=zeta_N(j-1-ki_est_i,N);%已知分数多普勒的时候多普勒维度的扩展
                end
                [maxze1,ze1index]=max(ze);
                %画图
                % figure()
                % plot(1:N,abs(ze),'b');
                % hold on;grid on;
                % plot(1:N, abs(Y(i,:))/Xp,'r');
                % legend("zeta","接收数据");
                h_est_i=abs(max(Y(i,:)))/abs((maxze1))/Xp;
                h_est=[h_est,h_est_i];
                temp=r(i, max_indices(i));
                h_phi_est_i=temp/abs(temp);
                h_phi_est=[h_phi_est,h_phi_est_i];
            end
        end
        % disp(['不同路径的时延=',num2str(li.')]);
        % disp(['不同路径的多普勒=',num2str(ki.')]);
        % disp(['不同路径的信道增益=',num2str(h_p_round)]);
        % disp(['不同路径的信道增益相位=',num2str(h_exp_round)]);
        % disp(['估计的时延抽头=',num2str(li_est)]);
        % disp(['估计的多普勒抽头=',num2str(ki_est)]);
        % disp(['估计的信道增益=',num2str(h_est)]);
        % disp(['估计的信道相位=',num2str(h_phi_est)]);



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
        num = norm(hw_est - hw, 'fro')^2;    % Σ_{l,k} |hw_est(l,k)-hw(l,k)|^2
        den = norm(hw, 'fro')^2;         % Σ_{l,k} |hw(l,k)|^2
        nmse(nums) = num / den;

        %% H_eff等效信道
        % y_vec = reshape(Y, M*N, 1);
        % dd_vec = reshape(dd, M*N, 1);
        % H_eff = zeros(M*N, M*N);
        % for l = 0:M-1       % 遍历矩阵的行
        %     for k = 0:N-1   % 遍历矩阵的列
        %         % 计算 Y(l+1, k+1) 在列主序向量 y_vec 中的索引
        %         row_idx = k*M + (l+1);
        %         for l_prime = 0:M-1
        %             for k_prime = 0:N-1
        %                 % 计算 dd(l_prime+1, k_prime+1) 在列主序向量 dd_vec 中的索引
        %                 col_idx = k_prime*M + (l_prime+1);
        %                 h_l = mod(l - l_prime, M);
        %                 h_k = mod(k - k_prime, N);
        %                 H_eff(row_idx, col_idx) = hw_est(h_l + 1, h_k + 1);
        %             end
        %         end
        %     end
        % end
        % y=H_eff*dd_vec;
        % yabs=reshape(abs(y),M,[]);
        % %% --- 导频干扰消除 (PIC) ---
        % 
        % % 步骤 1: 重建一个只包含导频的发送矩阵
        % % 创建一个全零矩阵，只在导频位置(1,1)放置导频值Xp
        % dd_pilot_only = zeros(M, N);
        % dd_pilot_only(1, 1) = Xp;
        % 
        % % 步骤 2: 估计导频部分对接收信号的贡献
        % % 通过将纯导频信号与我们估计的信道hw_est进行二维循环卷积
        % % 使用高效的FFT方法来完成这个卷积
        % Y_pilot_est = ifft2(fft2(dd_pilot_only) .* fft2(hw_est));
        % 
        % % 步骤 3: 从总接收信号中减去估计的导频影响
        % % Y_data_only 就是消除了导频干扰后，留给数据检测器处理的信号
        % Y_data_only = Y - Y_pilot_est;
        % 
        % 
        % %% MMSE 检测 (高效频域法)
        % H_lmmse = (H_eff' * H_eff + 1 * eye(size(H_eff, 2))) \ H_eff';
        % dd_est_matrix=H_lmmse*Y_data_only(:);
        % dd_est_matrix= reshape(dd_est_matrix, M, []);
        % 
        % % --- 步骤 5: 提取数据符号并解调 ---
        % % 找到原始数据符号的位置 (排除保护带的0和导频Xp)
        % data_indices = find(dd ~= 0 & dd ~= Xp);
        % 
        % % 提取对应位置的恢复符号
        % dd_est_vec = dd_est_matrix(data_indices);
        % dd_vec=dd(data_indices);
        % t1=abs(dd_est_vec);t2=abs(dd_vec);
        % 
        % dd_est_vec=dd_est_vec./sqrt(SNR/ans2);
        % % 提取原始发送的比特 (用于比较)
        % bit_orig = bit_all(data_indices);
        % 
        % % 解调恢复的符号，得到估计的比特
        % bit_est = qamdemod(dd_est_vec, modu, 'UnitAveragePower', true);
        % 
        % % --- 步骤 6: 计算本次迭代的比特误码率 (BER) ---
        % [num_errors, ber_iter(nums)] = biterr(bit_orig, bit_est);
        % 
        % [num_errors_this_frame, ~] = biterr(bit_orig(:), bit_est(:));
        % num_bits_this_frame = length(bit_orig(:));
        % total_bit_errors = total_bit_errors + num_errors_this_frame;
        % total_bits_transmitted = total_bits_transmitted + num_bits_this_frame;
        % if mod(nums, 50) == 0 % 每50次迭代打印一次状态
        %     fprintf('SNR: %.1f dB | Iter: %d | Total Errors: %d | Current BER: %e\n', ...
        %         SNR_dB(i_snr), nums, total_bit_errors, total_bit_errors/total_bits_transmitted);
        % end
    end
    chanel_err(i_snr)=mean(nmse);
    % ber(i_snr) = total_bit_errors / total_bits_transmitted;;
end
figure();% plot(SNR_dB,chanel_err);
semilogy(SNR_dB,chanel_err);
grid on;
%
% figure;
% semilogy(SNR_dB, ber, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
% xlabel('信噪比 (dB)');ylabel('比特误码率 (BER)');title('MMSE检测后系统BER性能曲线');legend('MMSE-BER');

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


