clc;clear;close all;

M=16;
N=16;
l_max=2;k_max=2;
P=2;%路径数
cen_modu = 2;       
guard_modu= 2;
bits_per_qam_symbol = log2(cen_modu);
bits_per_guard_symbol=log2(guard_modu);
g = 4; % IM分组参数
guard_data_power_factor =0.25;

%% 速率=位置信息+调制阶数
mo=mod((N-4*k_max-1),g);
q1 = floor((N-4*k_max-1)/ g);
%一行的速率与绿色区间激活个数
bit_row=(log2(guard_modu)+mo)+q1*(log2(guard_modu)+log2(nchoosek(g,1)));
bit_row_num=(1+q1)*(2*l_max+1);%guard个数
cen_num=(N/g)*(N-(1+2*l_max));
bit=(log2(cen_modu)+log2(nchoosek(g,1)))*(cen_num)+(2*l_max+1)*bit_row;

%% 信道检测
SNR_dB = 0:5:25;
iter=[2000,2000,2000,2000,2000,2000];
% SNR_dB=20;nmse_iter_per_snr=[2000];

%% 符号检测SNR
% SNR_dB=0:5:15;
% iter=[4e3,4e3,4e4,4e5];%bit数量

ber = zeros(1, length(SNR_dB));chanel_err=zeros(1,length(SNR_dB));
nmse_results = zeros(1, length(SNR_dB));

for i_snr=1:length(SNR_dB)
    SNR = 10^(SNR_dB(i_snr)/10);  
    total_bits_transmitted = 0;%传输bit
    total_bit_errors = 0;%错误bit
    nmse_accumulator = 0;%nmse累计
    nums = 0; % 当前迭代次数计数器

    for nums=1:iter(i_snr)%每次都重新生成信道
        nums
    % while  (nums*bit )<iter(i_snr)
    %     nums = nums + 1; % 迭代次数加1

        % 生成随机的延迟和分数多普勒值
        li = randperm(l_max, P)';
        k = randi([-k_max,k_max],P,1);kv=randi([-5,5],P,1)*0.1;ki=k+kv;
        % EVA信道，时延为0的路径是增益最大的，以后依次递减
        h_p_db=[1.5,1.4];
        h_p = 10.^((-1) * h_p_db/10)*10;%假设第一个路径增益为10
        h_p_round  = round(h_p, 4);
        h_exp=exp( 1j*2 * pi * rand(1, P));%生成路径相位
        h_exp_round = round(h_exp, 4);
        noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声
        noise=reshape(noise,M,N);
       

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
        

        %% 功率分配

        data_zone_mask = false(M,N);
        data_rows = (l_max+2):(M-l_max); % 数据行号
        data_zone_mask(data_rows, :) = true;

        % 找到在这个区域中被IM激活的符号:中央data
        active_data_indices = find(dd ~= 0 & data_zone_mask);
        num_central_symbols = length(active_data_indices);
        num_guard_symbols = ans3-num_central_symbols;
        % 解方程
        % 总数据功率 = 中央符号数*中央单位功率(默认为1) + 卫兵符号数*卫兵单位功率（1缩小guard_data_power_factor倍）
        denominator = num_central_symbols + num_guard_symbols * guard_data_power_factor;
        if denominator > 0
            power_per_central_symbol = SNR / denominator;
        else
            power_per_central_symbol = 0; % 防止除以0
        end
        power_per_guard_symbol = power_per_central_symbol * guard_data_power_factor;

        % 3.4 计算幅度缩放因子 (幅度是功率的平方根)
        scale_central = sqrt(power_per_central_symbol);
        scale_guard   = sqrt(power_per_guard_symbol);
       

        % 3.5 将不同的缩放因子应用到dd矩阵的不同区域
        % 注意：此时dd中的激活符号幅度还是1
        dd(data_rows,:) = dd(data_rows,:) * scale_central;
        guard_data_indices=[1:(l_max+1), (M-l_max+1):M];
        dd(guard_data_indices,:)   = dd(guard_data_indices,:) * scale_guard;

        % --- C: 放置导频并进行功率分配 ---
        Xp=10;    dd(1,1)=Xp;

        % 分数多普勒域下，DD域等效信道代码,发送数据和接受数据都是dd域
        hw=zeros(M,N);Y=zeros(M,N);
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
        Y=Y+noise;
        n=abs(noise);
        YABS=abs(Y);

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
                for i=1:P
                    theta=exp(1j*2*pi*ki_est(i)*(li_est(i))/(M*N));
                    delta_term = (l == li_est(i)); % 如果 l 等于 li(i)，则 delta 为 1,l是从0开始计算
                    hw_est(l+1,k+1)=hw_est(l+1,k+1)+h_est(i)*h_phi_est(i)*delta_term*zeta_N(k-ki_est(i),N).*theta;
                end
            end
        end
       
            % 如果变量 myVar 存在，则执行这里的代码
            num = norm(hw_est - hw, 'fro')^2;    % Σ_{l,k} |hw_est(l,k)-hw(l,k)|^2
            den = norm(hw, 'fro')^2;         % Σ_{l,k} |hw(l,k)|^2
            nmse(nums) = num / den;  
       

        
  % %% H_eff等效信道
  %       y_vec = reshape(Y, M*N, 1);
  %       dd_vec = reshape(dd, M*N, 1);
  %       H_eff = zeros(M*N, M*N);
  %       for l = 0:M-1       % 遍历矩阵的行
  %           for k = 0:N-1   % 遍历矩阵的列
  %               % 计算 Y(l+1, k+1) 在列主序向量 y_vec 中的索引
  %               row_idx = k*M + (l+1);
  %               for l_prime = 0:M-1
  %                   for k_prime = 0:N-1
  %                       % 计算 dd(l_prime+1, k_prime+1) 在列主序向量 dd_vec 中的索引
  %                       col_idx = k_prime*M + (l_prime+1);
  %                       h_l = mod(l - l_prime, M);
  %                       h_k = mod(k - k_prime, N);
  %                       H_eff(row_idx, col_idx) = hw_est(h_l + 1, h_k + 1);
  %                   end
  %               end
  %           end
  %       end
  %       y=H_eff*dd_vec;
  %       yabs=reshape(abs(y),M,[]);
  % 
  %       % --- 导频干扰消除 (PIC) ---
  % 
  %       % 步骤 1: 重建一个只包含导频的发送矩阵        % 创建一个全零矩阵，只在导频位置(1,1)放置导频值Xp
  %       dd_pilot_only = zeros(M, N);
  %       dd_pilot_only(1, 1) = Xp;
  %       % 步骤 2: 估计导频部分对接收信号的贡献       % 通过将纯导频信号与我们估计的信道hw_est进行二维循环卷积
  %       % 使用高效的FFT方法来完成这个卷积
  %       Y_pilot_est = ifft2(fft2(dd_pilot_only) .* fft2(hw_est));
  %       % 步骤 3: 从总接收信号中减去估计的导频影响        % Y_data_only 就是消除了导频干扰后，留给数据检测器处理的信号
  %       Y_data_only = Y - Y_pilot_est;
  % 
  %       % MMSE 检测 (高效频域法)
  %       H_lmmse = (H_eff' * H_eff + 1 * eye(size(H_eff, 2))) \ H_eff';
  %       dd_est_matrix=H_lmmse*Y_data_only(:);
  %       dd_est_matrix= reshape(dd_est_matrix, M, []);
  % 
  %       estimated_bits_stream=[];
  % 
  %       % ... 此前的代码已经计算出了估计矩阵 dd_est_matrix ...
  % 
  %       % 初始化一个空数组，用于存储在接收端检测出的比特流
  %       estimated_bits_stream = [];
  % 
  %       % --- A: 对中心数据区域进行检测 ---
  %       % 我们遍历发射端已知的网格结构，来检测传输的比特。
  %       for row = data_rows
  %           for j = 1:num_groups_data
  %               % 定义当前索引调制(IM)组所占的列
  %               start_col = (j - 1) * g + 1;
  %               group_cols = start_col : start_col + g - 1;
  %               % 从原始发送矩阵 dd 中提取当前组的符号，以找到真实激活位置
  %               original_group_symbols = dd(row, group_cols);
  %               % 【关键修改】: 直接使用真实位置，不再搜索最大值                        
  %               [~, true_local_idx] = max(abs(original_group_symbols));      
  %               im_bits_est = de2bi(true_local_idx - 1, bits_per_IM_data, 'left-msb');             
  %               group_est_symbols = dd_est_matrix(row, group_cols) / scale_central;
  %               symbol_to_demod = group_est_symbols(true_local_idx);
  %               qam_bits_est = qamdemod(symbol_to_demod, cen_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
  %               estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
  %           end
  %       end
  % 
  %       % --- B: 对保护间隔区域进行检测 ---
  %       for row = guard_rows
  %           i_col = 1;
  %           while i_col <= num_guard_cols
  %               current_group_size = min(g, num_guard_cols - i_col + 1);
  %               if current_group_size > 1
  %                   global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
  %                   % 【关键修改】: 直接从原始 dd 矩阵中找到真实激活位置
  %                   original_group_symbols = dd(row, global_col_indices);
  %                   [~, true_local_idx] = max(abs(original_group_symbols));
  %                   % 1. 索引检测 (理想化)
  %                   if current_group_size == g
  %                       num_im_bits = floor(log2(nchoosek(g, 1)));
  %                       im_bits_est = de2bi(true_local_idx - 1, num_im_bits, 'left-msb');
  %                   else
  %                       im_bits_est = zeros(1, current_group_size);
  %                       im_bits_est(true_local_idx) = 1;
  %                   end
  %                   % 2. 符号检测             
  %                   group_est_symbols = dd_est_matrix(row, global_col_indices) / scale_guard;
  %                   symbol_to_demod = group_est_symbols(true_local_idx);
  %                   qam_bits_est = qamdemod(symbol_to_demod, guard_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
  %                   estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
  %               end
  %               i_col = i_col + current_group_size;
  %           end
  %       end
  % 
  %       % % --- A: 处理中央数据区域 ---
  %       % % 检测中央数据区域
  %       % for row = data_rows
  %       %     for j = 1:num_groups_data
  %       %         start_col = (j - 1) * g + 1;               
  %       %         group_est = dd_est_matrix(row, start_col : start_col + g - 1) / scale_central;                
  %       %         [~, detected_local_idx] = max(abs(group_est));
  %       %         detected_symbol = group_est(detected_local_idx);
  %       %         im_bits_est = de2bi(detected_local_idx - 1, bits_per_IM_data, 'left-msb');
  %       %         qam_bits_est = qamdemod(detected_symbol, cen_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
  %       %         estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
  %       %     end
  %       % end
  %       % len1=length(estimated_bits_stream);
  %       % [~,err1]=biterr(original_bits_stream(1:len1), estimated_bits_stream);
  %       % 
  %       % % --- B: 处理保护间隔区域 ---
  %       % for row = guard_rows
  %       %     i_col = 1;
  %       %     while i_col <= num_guard_cols
  %       %         current_group_size = min(g, num_guard_cols - i_col + 1);
  %       %         if current_group_size  > 1 %与发射端设计呼应
  %       %             if current_group_size == g
  %       %                 global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
  %       %                 group_est = dd_est_matrix(row, global_col_indices) / scale_guard;
  %       %                 [~, detected_local_idx] = max(abs(group_est));
  %       %                 detected_symbol = group_est(detected_local_idx);
  %       %                 im_bits_est = de2bi(detected_local_idx - 1, bits_per_IM_data, 'left-msb');
  %       %                 qam_bits_est = qamdemod(detected_symbol, guard_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
  %       %                 estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
  %       %             else
  %       %                 global_col_indices = guard_cols_range(i_col : i_col + current_group_size - 1);
  %       %                 group_est = dd_est_matrix(row, global_col_indices) / scale_guard;
  %       %                 [~, detected_local_idx] = max(abs(group_est));
  %       %                 detected_symbol = group_est(detected_local_idx);
  %       %                 im_bits_est = zeros(1, current_group_size);
  %       %                 im_bits_est(detected_local_idx) = 1;
  %       %                 qam_bits_est = qamdemod(detected_symbol, guard_modu, 'UnitAveragePower', true, 'OutputType', 'bit');
  %       %                 estimated_bits_stream = [estimated_bits_stream, im_bits_est, qam_bits_est(:).'];
  %       %             end
  %       %         end
  %       %         i_col = i_col + current_group_size;
  %       %     end
  %       % end       
  %       % [~,err2]=biterr(original_bits_stream(len1+1:length(estimated_bits_stream)), estimated_bits_stream(len1+1:length(estimated_bits_stream)));
  %       %% 调试
  %       % A=reshape(original_bits_stream(len1+1:length(estimated_bits_stream)),bit_row,[]);
  %       % B=reshape(estimated_bits_stream(len1+1:length(estimated_bits_stream)),bit_row,[]);
  %       % IM1=q1*(log2(guard_modu)+log2(nchoosek(g,1)));
  %       % 
  %       % [num1,~]=biterr(A(1:IM1,:),B(1:IM1,:));%正常调制
  %       % [num2,~]=biterr(A(IM1+1:bit_row,:),B(IM1+1:bit_row,:));%独热编码
  %       % [~, ber_iter(nums)] = biterr(original_bits_stream, estimated_bits_stream);
  % 
  %       [num_errors_this_frame, ~] = biterr(original_bits_stream, estimated_bits_stream);
  %       num_bits_this_frame = length(original_bits_stream);
  %       total_bit_errors = total_bit_errors + num_errors_this_frame;
  %       total_bits_transmitted = total_bits_transmitted + num_bits_this_frame;
  % 
  % 
  % 
  %       if mod(nums, 50) == 0 % 每50次迭代打印一次状态
  %           fprintf('SNR: %.1f dB | Iter: %d | progress: %.2f%% | Current BER: %e\n', ...
  %                   SNR_dB(i_snr), nums, roundn(total_bits_transmitted/iter(i_snr),-2)*100, total_bit_errors/total_bits_transmitted);
  %       end

    end
    chanel_err(i_snr)=mean(nmse);
    % ber(i_snr) = total_bit_errors / total_bits_transmitted;

end

figure();
semilogy(SNR_dB,chanel_err);
grid on;hold on;
% figure();
% semilogy(SNR_dB, ber, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
% grid on;
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

