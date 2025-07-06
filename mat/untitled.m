clc;clear;close all;
M=16;
N=16;
hw=zeros(M,N);
P=1;
dd=zeros(M,N);
dd(1,:)=10;

Y=zeros(M,N);
li=0;ki=0.1;
h_p_round=1;

% g=4;%IM调制，每行分成几组
% dd(1,:)=10;
% for j = 1:N/g
%         % 计算当前组的起始列和结束列
%         start_col = (j-1) * (N/g) + 1;
%         end_col = j * (N/g);
%         % 每组只有一个激活
%         zero_indices = start_col:end_col;
%         num_zeros = 3;
%         zero_positions = randperm(length(zero_indices), num_zeros);
%         % 将随机选择的部分置为 0
%         dd(1, zero_indices(zero_positions)) = 0;
% end


for l=0:M-1
    for k=0:N-1
        for i=1:P
            theta=exp(1j*2*pi*ki(i)*(li(i))/(M*N));
            delta_term = (l == li(i)); % 如果 l 等于 li(i)，则 delta 为 1
            hw(l+1,k+1)=hw(l+1,k+1)+h_p_round(i)*delta_term*zeta_N(k-ki(i),N).*theta;
        end
    end
end

       y_vec = reshape(Y, M*N, 1);
       dd_vec = reshape(dd, M*N, []);
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
        noise=reshape(noise,M,N);
        y=H_eff*dd_vec;
        Y=reshape(y,M,[])+noise;
        H_lmmse = (H_eff' * H_eff + 1 * eye(size(H_eff, 2))) \ H_eff';
        dd_est_matrix=H_lmmse*Y(:);
        dd_est_matrix= reshape(dd_est_matrix, M, []);
        b=abs(dd_est_matrix);

        figure();
        plot(1:16,abs(dd(1,:)),'-o');
        hold on;
        grid on;
        plot(1:16,b(1,:),'-^');
        legend('发射信号','均衡的接收信号');ylabel('幅值');
        xlabel('满铺');

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


