clc;clear;close all;

M=16;
N=16;
fc=77e9;%GHz
deltaf=195e3;%KHz
T=1/deltaf;
tau_resolution=1/(M*deltaf);
doppler_resolution=1/(N*T);
l_max=2;
k_max=2;
P=2;%路径数
iter=100;
nmse=zeros(1,iter);
doppler_err=zeros(1,iter);delay_err=zeros(1,iter);h_err=zeros(1,iter);h_phi_err=zeros(1,iter);
for nums=1:iter;
    nums
li = randperm(l_max, P)';
k = randi([-k_max,k_max],P,1);kv=randi([-5,5],P,1)*0.1;ki=k+kv;
for j=1:P
    if(ki<0) 
        ki=-fix(ki*10)/10;
    else
        ki=fix(ki*10)/10;
    end
end


% h设置为正实数，同时设置h的相位
% h_EVA_dB=[1.5,1.4,3.6,0.6,9.1,7,12,16.9];%EVA信道参数
h_exp=exp( 1j*2 * pi * rand(1, P));%生成增益相位
h_p_db=[1.5,1.4];
h_p = 10.^((-1) * h_p_db/10)*10;%假设第一个路径增益为10
noise=sqrt(1/2)*(randn(M*N,1)+1i*randn(M*N,1));%均值为0方差为1的高斯噪声

% disp(['不同路径的时延=',num2str(li.')]);
% disp(['不同路径的多普勒=',num2str(ki.')]);
% disp(['不同路径的信道增益=',num2str(h_p)]);
% disp(['不同路径的信道增益相位=',num2str(h_exp)]);


%TEST DATA
% li=[1];ki=[0];
%% 生成bit流
% modu=4;
% bit_all = randi(modu,M,N)-ones(M,N); %DD域
% x = qammod(bit_all,modu,'UnitAveragePower',true);%调制后的信号的平均功率为 1，总能量/个数，能量为模长平方
% dd=reshape(x,M,N);
% g=4;%IM调制，每行分成几组
% for i = 1:M
%     for j = 1:g
%         % 计算当前组的起始列和结束列
%         start_col = (j-1) * (N/g) + 1;
%         end_col = j * (N/g);
%         % 每组只有一个激活
%         zero_indices = start_col:end_col;
%         num_zeros = 3;
%         zero_positions = randperm(length(zero_indices), num_zeros);
%         % 将随机选择的部分置为 0
%         dd(i, zero_indices(zero_positions)) = 0;
%     end
% end
% 导频位置将导频放在1,1左上角
dd=zeros(M,N);
SNR_PILOT_dB=20;
Xp=10^(SNR_PILOT_dB/10);
% for i=1:l_max+1
%     for j=1:2*k_max+1
%         dd(i,j)=0;
%     end
% end
% for i=1:l_max+1
%     for j=N:-1:N-2*k_max
%         dd(i,j)=0;
%     end
% end
% for i=M-l_max+1:1:M
%     for j=1:2*k_max+1
%         dd(i,j)=0;
%     end
% end
% for i=M-l_max+1:M
%     for j=N:-1:N-2*k_max
%         dd(i,j)=0;
%     end
% end
dd(1,1)=Xp;
hw=zeros(M,N);
Y=zeros(M,N);

%% 分数多普勒域下，DD域等效信道代码,发送数据和接受数据都是dd域
for l=0:M-1
    for k=0:N-1
        for i=1:P
            theta=exp(1j*2*pi*ki(i)*(l-li(i))/(M*N));
            delta_term = (l == li(i)); % 如果 l 等于 li(i)，则 delta 为 1
            hw(l+1,k+1)=hw(l+1,k+1)+h_p(i)*h_exp(i)*delta_term*zeta_N(k-ki(i),N).*theta;
        end
    end
end
magnitude = abs(hw);
% figure();
% plot(0:M-1,magnitude(2,:),'bo');hold on;grid on;

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

% front_col= Y(:,1:2*k_max+1);
% back_col = Y(:,N-2*k_max+1:N);
% result0 = [back_col, front_col];%将保护间隔位置进行合并
% mag_result=abs(result);

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
    if(max_values(i)>3)%存在时延
        li_est=[li_est,i-1];
        ki_est_i=(-k_max-0.5)+step*(max_indices(i)-1);%确定时延以后，得到对应时延的多普勒数值
        ki_est=[ki_est,ki_est_i]; 
        ze=zeros(1,N);
        for j=1:N
        ze(j)=zeta_N(j-1-ki_est_i,N);%已知分数多普勒的时候多普勒维度的扩展      
        end
        [maxze1,ze1index]=max(ze);
        % figure()
        % plot(1:N,abs(ze),'b');
        % hold on;grid on;
        % plot(1:N, abs(Y(i,:))/Xp,'r');
        % legend("zeta","接收数据");
        h_est_i=abs(max(Y(i,:)))/abs((maxze1))/Xp;
        h_est=[h_est,h_est_i];
        h_phi_est_i=max(Y(i,:))/(maxze1)/Xp/h_est_i  ; 
        h_phi_est=[h_phi_est,h_phi_est_i];
    end
end
% 
% disp(['估计的时延抽头=',num2str(li_est)]);
% disp(['估计的多普勒抽头=',num2str(ki_est)]);
% disp(['估计的信道增益=',num2str(h_est)]);
% disp(['估计的信道相位=',num2str(h_phi_est)]);

delay_err(nums)=sum(sort(li)-sort(li_est.'));
% if()<1e-4)
doppler_err(nums)=sum(sort(ki)-sort(ki_est.'));%转置.'
h_err(nums)=sum(sort(h_p)-sort(h_est));
h_phi_err(nums)=sum(sort(h_exp)-sort(h_phi_est));

% hw_est=zeros(M,N);
% for l=0:M-1
%     for k=0:N-1
%         for i=1:P
%             theta=exp(1j*2*pi*ki_est(i)*(l-li_est(i))/(M*N));
%             delta_term = (l == li_est(i)); % 如果 l 等于 li(i)，则 delta 为 1
%             hw_est(l+1,k+1)=hw_est(l+1,k+1)+h_est(i)*h_phi_est(i)*delta_term*zeta_N(k-ki_est(i),N).*theta;
%         end
%     end
% end
% h_error=hw_est-hw;
% nmse(nums)=norm(h_error,2)/norm(hw,2);
end
sum(delay_err)
sum(doppler_err)
sum(h_err)
abs(sum(h_phi_err))
% figure();plot(1:iter,abs(h_err));figure();plot(1:iter,abs(h_phi_err));
% sum(nmse)/iter


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

%% 绘制图形
% k_values = linspace(0, N, 1000);  % 生成从0到N的1000个点
% zeta_values = arrayfun(@(k) zeta_N(k-0.4, N), k_values);  % 计算 zeta_N(k)figure;
% % figure(2);
% plot(k_values, abs(zeta_values), 'Color', [1, 0.647, 0]);  % 使用橙色绘制曲线
% grid on; result(l,:)
