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
%路径数
P=1;
li = randi(l_max,P,1);
k = randi([-k_max,k_max],P,1);
kv=randi([-5,5],P,1)*0.1;
ki=k+kv;

li=[1];
ki=[-1.8];
% h=[1.5,3.6];
h=[1,1];
%生成bit流
modu=4;
bit_all = randi(modu,M,N)-ones(M,N); %DD域
x = qammod(bit_all,modu,'UnitAveragePower',true);%调制后的信号的平均功率为 1，总能量/个数，能量为模长平方
dd=reshape(x,M,N);
%索引调制
g=4;%每行分成几组
for i = 1:M
    for j = 1:g
        % 计算当前组的起始列和结束列
        start_col = (j-1) * (N/g) + 1;
        end_col = j * (N/g);
        % 每组只有一个激活
        zero_indices = start_col:end_col;
        num_zeros = 3;
        zero_positions = randperm(length(zero_indices), num_zeros);
        % 将随机选择的部分置为 0
        dd(i, zero_indices(zero_positions)) = 0;
    end
end

dd=zeros(M,N);
%poilt位置位于0,0
Xp=1;
for i=1:l_max+1
    for j=1:2*k_max+1
        dd(i,j)=0;
    end
end
for i=1:l_max+1
    for j=N:-1:N-2*k_max
        dd(i,j)=0;
    end
end
for i=M-l_max+1:1:M
    for j=1:2*k_max+1
        dd(i,j)=0;
    end
end
for i=M-l_max+1:M
    for j=N:-1:N-2*k_max
        dd(i,j)=0;
    end
end
dd(1,1)=Xp;

hw=zeros(M,N);
Y=zeros(M,N);


%% DD域等效信道代码,发送数据和接受数据都是dd域
for l=0:M-1
    for k=0:N-1
        for i=1:P
            theta=exp(1j*2*pi*ki(i)*(l-li(i))/(M*N));
            delta_term = (l == li(i)); % 如果 l 等于 li(i)，则 delta 为 1
            hw(l+1,k+1)=hw(l+1,k+1)+h(i)*delta_term*zeta_N(k-ki(i),N).*theta;
        end
    end
end
figure(1);
magnitude = abs(hw);
plot(0:M-1,magnitude(2,:),'bo');hold on;grid on;


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
[a,b]=max(Y(2,:));
y2abs=abs(Y(2,:));
plot(0:M-1,y2abs,'r*');hold on;

front_col= Y(:,1:2*k_max+1);
back_col = Y(:,N-2*k_max+1:N);
result = [back_col, front_col];
mag_result=abs(result);
step=0.1;

yd=zeros(1,4*k_max+1);
colnum=((2*k_max+1)/(step))+1;
r=zeros(l_max+1,colnum);
for l=1:l_max+1
    j=1;
    yd=result(l,:);
    for stepval=-k_max-0.5:step:k_max+0.5
        for k=1:4*k_max+1
            r(l,j)=r(l,j)+yd(k)*conj(zeta_N(k-stepval-2*k_max-1,N));           
        end
         j=j+1;
    end
end
figure(2);
ydd=zeros(4*k_max+1,1);
ze=zeros(4*k_max+1,1);
it=1:4*k_max+1;
for i=1:4*k_max+1
    ydd(i)=result(2,i);
    ze(i)=zeta_N(i-ki-2*k_max-1,N);
end
plot(it,abs(ze),'b');hold on;grid on;
plot(it,abs(ydd),'r');

rabs=abs(r);
[max_values, max_indices] = max(rabs, [], 2);
ki_est(i)=-k_max-0.5+step*(max_indices(2)-1);




%% 绘制图形
% k_values = linspace(0, N, 1000);  % 生成从0到N的1000个点
% zeta_values = arrayfun(@(k) zeta_N(k-0.4, N), k_values);  % 计算 zeta_N(k)figure;
% % figure(2);
% plot(k_values, abs(zeta_values), 'Color', [1, 0.647, 0]);  % 使用橙色绘制曲线
% grid on; result(l,:)

function output = zeta_N(k, N)
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


