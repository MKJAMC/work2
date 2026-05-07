clc; clear; close all;
% 需要划分为4000的训练集，400的验证集，200的测试集
% 确保这里的文件名是你刚刚生成了三通道数据的那个文件名
load('data4600.mat'); 

%% ================= 维度检查 (新增) =================
% 检查加载的数据是否包含 3 个通道
[M, N, C, Batch] = size(Dataset_Y);
fprintf('成功加载数据，维度为: [%d, %d, %d, %d]\n', M, N, C, Batch);

if C ~= 3
    warning('注意：当前数据的通道数不是 3 (实部+虚部+模长)，而是 %d。请确认上一步生成是否正确。', C);
else
    fprintf('确认：检测到 3 通道数据 (实部, 虚部, 模长)。\n');
end

%% ================= 第三部分：数据集顺序划分与保存 =================
fprintf('3. 正在按顺序进行数据集划分...\n');

% --- 1. 定义索引范围 ---
% 训练集: 1 ~ 4000
% 验证集: 4001 ~ 4400
% 测试集: 4401 ~ 4600
train_range = 1:4000;
val_range   = 4001:4400;
test_range  = 4401:4600;

% --- 2. 划分输入数据 X (维度: M x N x 3 x Batch) ---
% 注意：这里的第三维用 : 表示全选，因此它会自动适配 3 通道
X_train = Dataset_Y(:, :, :, train_range);
X_val   = Dataset_Y(:, :, :, val_range);
X_test  = Dataset_Y(:, :, :, test_range);

% --- 3. 准备并划分标签数据 Y ---
% 合并时延和多普勒索引
% Target_Delay_Index: [4600, 3]
% Target_Doppler_Index: [4600, 3]
% cat(3, ...) 后的维度: [4600, 3, 2] (Batch, Targets, Labels)
Labels_All = cat(3, Target_Delay_Index, Target_Doppler_Index);

Y_train = Labels_All(train_range, :, :);
Y_val   = Labels_All(val_range, :, :);
Y_test  = Labels_All(test_range, :, :);

% --- 4. 提取并划分 SNR 标签 ---
% 确保生成代码中保存了 Target_SNR 变量
if exist('Target_SNR', 'var')
    SNR_train = Target_SNR(train_range);
    SNR_val   = Target_SNR(val_range);
    SNR_test  = Target_SNR(test_range);
    
    fprintf('划分完成（顺序截取）：\n');
    fprintf('  - 训练集: 1-4000 (SNR: %d-%d dB)\n', min(SNR_train), max(SNR_train));
else
    warning('未找到 Target_SNR 变量，将跳过 SNR 标签划分。');
    SNR_train = []; SNR_val = []; SNR_test = [];
end

% --- 5. 保存数据 (包含 SNR 标签) ---
% 保存为 -v7.3 格式以供 Python h5py 读取
save('OTFS_Dataset_Sequential.mat', ...
    'X_train', 'Y_train', 'SNR_train', ...
    'X_val',   'Y_val',   'SNR_val',   ...
    'X_test',  'Y_test',  'SNR_test',  '-v7.3');

fprintf('数据集已保存至 OTFS_Dataset_Sequential.mat\n');
fprintf('  X_train 最终维度: %s\n', mat2str(size(X_train)));
fprintf('  (Python读取时将变为 Batch x 3 x N x M)\n');