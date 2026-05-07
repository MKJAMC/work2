clc;
clear;
close all;

%% 1. 参数设置与数据生成
% 为了方便观察，我们使用一个较小的矩阵尺寸
M = 4; % 矩阵的行数
N = 3; % 矩阵的列数

% 创建一个 M x N 的随机复数矩阵 Y
Y = rand(M, N) + 1i * rand(M, N);

fprintf('将要验证的核心关系是: \n');
fprintf('vec(fft2(Y))  <==>  (F_N ⊗ F_M) * vec(Y) \n\n');

%% 2. 方法一：使用高效算法 fft2 (计算等式左手边 LHS)
% 这是在实际编程中我们使用的方法
Y_freq = fft2(Y);
lhs = Y_freq(:); % 将 M x N 的结果矩阵按列拉直成 MN x 1 的向量

fprintf('--- 方法一 (fft2 算法) 的结果 (部分显示) ---\n');
disp(lhs(1:5)); % 只显示前5个元素

%% 3. 方法二：使用线性代数定义 (计算等式右手边 RHS)
% 这是在理论推导中我们使用的方法

% a) 创建 M x M 的一维DFT矩阵 F_M
% dftmtx 函数可以直接生成标准的DFT矩阵
F_M = dftmtx(M);

% b) 创建 N x N 的一维DFT矩阵 F_N
F_N = dftmtx(N);

% c) 使用克罗内克积 (Kronecker Product) 构建 MN x MN 的二维DFT矩阵 F
% 注意 kron(A, B) 对应 A ⊗ B
F = kron(F_N, F_M);

% d) 将原始矩阵 Y 按列拉直成向量 y_vec
y_vec = Y(:);

% e) 执行矩阵-向量乘法
rhs = F * y_vec;

fprintf('--- 方法二 (线性代数定义) 的结果 (部分显示) ---\n');
disp(rhs(1:5)); % 只显示前5个元素

%% 4. 验证两者结果
% 计算两个结果向量之间的差值
difference = lhs - rhs;

% 计算差值向量的范数(长度)。如果两者相等，范数应该是一个非常接近0的数
% (由于计算机浮点精度，通常不会是绝对的0)
error_norm = norm(difference);

fprintf('--- 验证结果 ---\n');
fprintf('两种方法计算结果的差值范数: %e\n', error_norm);

% 设置一个很小的阈值来判断是否相等
if error_norm < 1e-10
    fprintf('\n验证成功！🎉\n');
    fprintf('这证明了 fft2(Y) 的拉直结果与 F * Y(:) 的数学定义完全等价。\n');
else
    fprintf('\n验证失败。\n');
end