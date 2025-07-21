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
k_max0=doppler_max*N*(1/delta_f);
% EVA信道，时延为0的路径是增益最大的，以后依次递减
h_p_db=[1.5,1.4,3.6,0.6,9.1,7,12,16.9];
relative_power_linear = 10.^((-1) * h_p_db/10);%db转换为功率
total_power = sum(relative_power_linear);
h_p = relative_power_linear / total_power;%每一个增益所占据的百分比
h_exp=exp( 1j*2 * pi * rand(1, length(h_p_db)));
l_max=20;k_max=2;
P=length(delay);%路径数
cen_modu =2;
guard_modu= 2;
bits_per_qam_symbol = log2(cen_modu);
bits_per_guard_symbol=log2(guard_modu);
g = 4; % IM分组参数
guard_data_power_factor =0.25;
%% 速率=位置信息+调制阶数
mo=mod((N-4*k_max-1),g);
q1 = floor((N-4*k_max-1)/ g);
%一行的速率与绿色区间激活个数
bit_row=(log2(guard_modu)+mo)+q1*(log2(guard_modu)+log2(nchoosek(g,1)));%guard一行bit数
bit_row_num=(1+q1)*(2*l_max+1);%guard个数
cen_num=(N/g)*(M-(1+2*l_max));
ans3=cen_num+bit_row_num;
bit=(log2(cen_modu)+log2(nchoosek(g,1)))*(cen_num)+(2*l_max+1)*bit_row;

%% 信道检测
SNR = 10:2:16;
iter=[2000,2000,2000,2000];
Xp_snr=60;
power_persym=10.^(SNR/10);%每个符号的功率
ans2=(M-2*l_max-1)*N+(N-4*k_max-1)*(2*l_max+1);
power=power_persym*ans2;%ans2为基准，求得总功率，方案二的功率为1
nmse_results = zeros(1, length(SNR));

for i_snr=1:length(SNR)

    nmse_accumulator = 0;%nmse累计
    nums = 0; % 当前迭代次数计数器
    %% 功率分配
   data_zone_mask = false(M,N);
    data_rows = (l_max+2):(M-l_max); % 数据行号
    data_zone_mask(data_rows, :) = true;
    % 找到在这个区域中被IM激活的符号:中央data
    num_central_symbols = cen_num;
    num_guard_symbols = ans3-num_central_symbols;
    % 解方程
    % 总数据功率 = 中央符号数*中央单位功率(默认为1) + 卫兵符号数*卫兵单位功率（1缩小guard_data_power_factor倍）
    denominator = num_central_symbols + num_guard_symbols * guard_data_power_factor;
    if denominator > 0
        power_per_central_symbol = power(i_snr) / denominator;
    else
        power_per_central_symbol = 0; % 防止除以0
    end
    power_per_guard_symbol = power_per_central_symbol * guard_data_power_factor;
    % 3.4 计算幅度缩放因子 (幅度是功率的平方根
    scale_central = sqrt(power_per_central_symbol);
    scale_guard   = sqrt(power_per_guard_symbol);

    for nums=1:iter(i_snr)%每次都重新生成信道
        nums

        %%  发射端：基于信息比特生成OTFS-IM信号
        dd = zeros(M, N);
        original_bits_stream = [];
        data_rows = (l_max + 2):(M - l_max);
        num_groups_data = N / g; % 计算分组数量
        bits_per_IM_data = log2(nchoosek(g, 1));%每组激活一个

        % --- A: 处理中央数据区域 (g = 组大小) ---
        for row = data_rows
            for j = 1:num_groups_data % 循环 N/g 次
                im_bits = randi([0 1], 1, bits_per_IM_data);%生成每一组的IM位置bit
                qam_bits = randi([0 1], 1, bits_per_qam_symbol);%生成一个qam符号的bit
                original_bits_stream = [original_bits_stream, im_bits, qam_bits];

                active_local_idx = bi2de(im_bits, 'left-msb') + 1;%将IM的bit转为位置
                symbol_int = bi2de(qam_bits, 'left-msb');
                qam_symbol = qammod(symbol_int, cen_modu, 'UnitAveragePower', true);
                start_col = (j - 1) * g + 1; % 每组的起始列
                dd(row, start_col + active_local_idx - 1) = qam_symbol;
            end
        end

        % --- B: 处理保护间隔区域 (使用独热编码, g = 组大小) ---
        guard_rows = [1:(l_max + 1), (M - l_max + 1):M];
        guard_cols_range = (2 * k_max + 2):(N - 2 * k_max);
        num_guard_cols = length(guard_cols_range);%需要索引调制的保护间隔部分

        % 遍历每一个指定的保护行
        for row = guard_rows
            i_col = 1; % 初始化保护列的起始计数器
            while i_col <= num_guard_cols     % 计算当前组的大小，最大不超过 g
                current_group_size = min(g, num_guard_cols - i_col + 1);
                if current_group_size > 1  %最后余数为1，则不存在索引调制
                    if current_group_size == g
                        num_im_bits = floor(log2(nchoosek(g, 1)));
                        im_bits = randi([0 1], 1, num_im_bits);
                        active_local_idx = bi2de(im_bits, 'left-msb') + 1;
                    else
                        % --- 情况2: 剩余组 (current_group_size < g)，使用独热编码,在当前剩余的组内随机选择一个激活位置
                        active_local_idx = randi(current_group_size);
                        % 生成独热编码比特向量
                        im_bits = zeros(1, current_group_size);
                        im_bits(active_local_idx) = 1;
                    end
                    qam_bits = randi([0 1], 1, bits_per_guard_symbol);
                    original_bits_stream = [original_bits_stream, im_bits, qam_bits];
                    % 将QAM比特调制成一个复数符号
                    symbol_int = bi2de(qam_bits, 'left-msb');
                    qam_symbol = qammod(symbol_int, guard_modu, 'UnitAveragePower', true);
                    % 获取当前组在dd的实际列索引
                    global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
                    % 将QAM符号放置在被选中的激活位置上
                    dd(row, global_col_indices(active_local_idx)) = qam_symbol;
                end
                % 移动到下一个组的起始位置
                i_col = i_col + current_group_size;
            end
        end

        ans3 = nnz(dd);%总激活个数

        % 3.5 将不同的缩放因子应用到dd矩阵的不同区域
        % 注意：此时dd中的激活符号幅度还是1
        dd(data_rows,:) = dd(data_rows,:) .* scale_central;
        guard_data_indices=[1:(l_max+1), (M-l_max+1):M];
        dd(guard_data_indices,:)   = dd(guard_data_indices,:) * scale_guard;

        % --- C: 放置导频并进行功率分配 ---
        Xp=sqrt(10^(Xp_snr/10));    dd(1,1)=Xp;

        % 分数多普勒域下，DD域等效信道代码,发送数据和接受数据都是dd域
        hw=zeros(M,N);Y=zeros(M,N);
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

        %% 添加噪声
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声
        noise=reshape(noise,M,N);
        Y=Y+noise;


        % 互相关矩阵
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
            if(max_values(i)>3)%判断时延是否存在的证据，三倍的标准差
                li_est=[li_est,i-1];
                ki_est_i=(-k_max-0.5)+step*(max_indices(i)-1);%确定时延以后，得到对应时延的多普勒数值
                ki_est=[ki_est,ki_est_i];
                ze=zeros(1,N);
                for j=1:N
                    ze(j)=zeta_N(j-1-ki_est_i,N);%已知分数多普勒的时候多普勒维度的扩展
                end
                [maxze1,ze1index]=max(ze);
                h_est_i=abs(max(Y(i,:)))/abs((maxze1))/Xp;
                h_est=[h_est,h_est_i];
                temp=r(i, max_indices(i));
                h_phi_est_i=temp/abs(temp);
                h_phi_est=[h_phi_est,h_phi_est_i];
            end
        end

        %% 估计的数值生成信道,NMSE
        hw_est=zeros(M,N);
        for l=0:M-1
            for k=0:N-1
                for i=1:length(li_est)
                    theta=exp(1j*2*pi*ki_est(i)*(li_est(i))/(M*N));
                    delta_term = (l == li_est(i)); % 如果 l 等于 li(i)，则 delta 为 1,l是从0开始计算
                    hw_est(l+1,k+1)=hw_est(l+1,k+1)+h_est(i)*h_phi_est(i)*delta_term*zeta_N(k-ki_est(i),N).*theta;
                end
            end
        end

        %% 估计的数值生成信道,NMSE
        hw_est=zeros(M,N);
        for l=0:M-1
            for k=0:N-1
                for i=1:length(li_est)
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
% grid on;hold on;


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

