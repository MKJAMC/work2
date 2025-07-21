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
theta0= -pi + 2*pi*rand(1, length(delay));
doppler=doppler_max*cos(theta0) ;
ki=doppler*N*(1/delta_f);
% % ki=round(ki,1);
k_max0=doppler_max*N*(1/delta_f);

% EVA信道，时延为0的路径是增益最大的，以后依次递减
h_p_db=[1.5,1.4,3.6,0.6,9.1,7,12,16.9];
relative_power_linear = 10.^((-1) * h_p_db/10);%db转换为功率
total_power = sum(relative_power_linear);
h_p = relative_power_linear / total_power;%每一个增益所占据的百分比
h_exp=exp( 1j*2 * pi * rand(1, length(h_p_db)));

l_max=20;k_max=2;
P=length(delay);%路径数
modu=2;%BPSK
% 方案一的数据个数
ans1=(M-2*l_max-1)*N;
% 方案一的速率
bit= (M-2*l_max-1)*N*log2(modu);
%方案二的数据个数,其MN,l_max,k_max不变则个数不变
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
%% 信道估计SNR
SNR = 10:2:16;
iter=[2000,2000,2000,2000];%迭代次数
Xp_snr=60;

power_persym = 10.^(SNR/10)*log2(modu); % 每个符号的能量，考虑到高阶调制
power=power_persym*ans2;%高阶调制后的全部符号能量
factor=power./ans2;%计算真正每个符号的能量

chanel_err=zeros(1,length(SNR));
% 方案一的数据个数
ans1=(M-2*l_max-1)*N;
% 方案一的速率
bit= (M-2*l_max-1)*N*log2(modu);
%方案二的数据个数,其MN,l_max,k_max不变则个数不变
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);


for i_snr=1:length(SNR)
    nmse=zeros(1,iter(i_snr));
    nums=0;
    for nums=1:iter(i_snr)%每次都重新生成信道
        nums

        %% 生成bit流
        bit_all = randi(modu,M,N)-ones(M,N); %DD域
        x = qammod(bit_all,modu,'UnitAveragePower',true);%调制后的信号的平均功率为 1，总能量/个数，能量为模长平方
        dd=reshape(x,M,N);
        dd=dd.*sqrt((factor(i_snr)));

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


        Xp=sqrt(10^(Xp_snr/10));    dd(1,1)=Xp;
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
        %% 噪声
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声
        noise=reshape(noise,M,N);
        Y=Y+noise;
        YABS=abs(Y);

        %% 互相关信道估计
        step=0.1;
        yd=zeros(1,size(Y,2));%待抽取出的每一行
        colnum=((2*k_max+1)/(step))+1;%step细化以后的总列数
        r=zeros(l_max+1,colnum);%相关值进行存放，行代表着带估计的时延位置，列代表着带估计的分数多普勒的位置
        for l=1:l_max+1
            j=1;
            yd=Y(l,:);
            for stepval=-k_max:step:k_max
                for k=1:size(Y,2)
                    r(l,j)=r(l,j)+yd(k)*conj(zeta_N(k-stepval-1,N));   %计算出不同step下的自相关值,-1是因为k是从1开始计算的
                end
                j=j+1;%不同的j代表不同的stepval
            end
        end
        % disp(['不同路径的时延=',num2str(li)]);
        % disp(['不同路径的多普勒=',num2str(ki)]);
        % disp(['不同路径的信道增益=',num2str(h_p)]);
        % disp(['不同路径的信道相位=',num2str(h_exp)]);

        li_est = []; h_est=[];h_phi_est=[];ki_est=[];
        rabs=abs(r);
        [max_values, max_indices] = max(rabs, [], 2);%选出了每一行的最大值

        for i=1:length(max_values)
            yd= Y(i,:);
            if(max(YABS(i,:))>2.8)%判断时延是否存在的证据
                li_est=[li_est,i-1];
                k_hat=(-k_max)+step*(max_indices(i)-1);%确定时延以后，得到对应时延的多普勒数值  

            %     peak_idx = max_indices(i); % 粗略搜索找到的峰值索引
            %     % 确定边界索引，并处理边界情况
            %     % 保证索引不越界 (例如峰值在第一个或最后一个点)
            %     left_idx = max(1, peak_idx - 2);
            %     right_idx = min(colnum, peak_idx + 2);
            %     % 使用峰值点和它两边的点作为 GSS 的初始边界
            %     b_l = (-k_max) + step * (left_idx - 1);
            %     b_u = (-k_max) + step * (right_idx - 1);                         
            %     max_iter=150;
            %     eta = (sqrt(5)-1)/2;
            %     % 在循环外计算初始的 b1 和 b2 以及它们的函数值
            %     b1 = b_u - eta * (b_u - b_l);
            %     b2 = b_l + eta * (b_u - b_l);
            %     f1_val = 0;
            %     for k=1:size(Y,2)
            %         f1_val = f1_val + yd(k)*conj(zeta_N(k-b1-1,N));
            %     end
            %     f1_val = abs(f1_val); % 取模
            % 
            %     f2_val = 0;
            %     for k=1:size(Y,2)
            %         f2_val = f2_val + yd(k)*conj(zeta_N(k-b2-1,N));
            %     end
            %     f2_val = abs(f2_val); % 取模
            % 
            %     for iter_gss = 1:max_iter
            %         if f1_val > f2_val
            %             b_u = b2;
            %             b2 = b1;
            %             f2_val = f1_val; % 无需重新计算 f2
            %             b1 = b_u - eta * (b_u - b_l);
            %             f1_val = 0;
            %             for k=1:size(Y,2)
            %                 f1_val = f1_val + yd(k)*conj(zeta_N(k-b1-1,N));
            %             end
            %             f1_val = abs(f1_val); % 只需要计算一个新的 f1
            %         else
            %             b_l = b1;
            %             b1 = b2;
            %             f1_val = f2_val; % 无需重新计算 f1
            %             b2 = b_l + eta * (b_u - b_l);
            %             f2_val = 0;
            %             for k=1:size(Y,2)
            %                 f2_val = f2_val + yd(k)*conj(zeta_N(k-b2-1,N));
            %             end
            %             f2_val = abs(f2_val); % 只需要计算一个新的 f2
            %         end
            %     end
            %     nu_opt = (b_l + b_u) / 2;
            % end
            % nu_opt-ki(1)
            % k_hat-ki(1)
           ki_est=[ki_est,k_hat];

            ze=zeros(1,N);
            for j=1:N
                ze(j)=zeta_N(j-1-k_hat,N);%已知分数多普勒的时候多普勒维度的扩展
            end
            [maxze1,ze1index]=max(ze);
            h_est_i=abs(max(Y(i,:)))/abs((maxze1))/Xp;
            h_est=[h_est,h_est_i];
            temp=r(i, max_indices(i));
            h_phi_est_i=temp/abs(temp);
            h_phi_est=[h_phi_est,h_phi_est_i];
        end

        % disp(['估计的时延抽头=',num2str(li_est)]);
        % disp(['估计的多普勒抽头=',num2str(ki_est)]);
        % disp(['估计的信道增益=',num2str(h_est)]);
        % disp(['估计的信道相位=',num2str(h_phi_est)]);
        end
        % delay_err(nums)=sum(sort(li)-sort(li_est));
        % doppler_err(nums)=sum(sort(ki) - sort(ki_est) );%转置.'
        % h_err(nums)=sum(sort(h_p)-sort(h_est));
        % h_phi_err(nums)=sum(sort(h_exp)-sort(h_phi_est));

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


    end
    chanel_err(i_snr)=mean(nmse);

end
figure();
% semilogy(SNR_dB,chanel_err);
grid on;xlabel("SNR");ylabel("NMSE");



