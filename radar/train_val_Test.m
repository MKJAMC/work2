clc; clear; close all;
% 需要划分为4000的训练集，400的验证集，200的测试集，分别保存到对应变量中
load('data4600.mat');

%% ================= 第三部分：数据集顺序划分与保存 =================
fprintf('3. 正在按顺序进行数据集划分...\n');

% --- 1. 定义索引范围 ---
% 训练集: 1 ~ 4000
% 验证集: 4001 ~ 4400
% 测试集: 4401 ~ 4600
train_range = 1:4000;
val_range   = 4001:4400;
test_range  = 4401:4600;

% --- 2. 划分输入数据 X (维度: M x N x 2 x Batch) ---
X_train = Dataset_Y(:, :, :, train_range);
X_val   = Dataset_Y(:, :, :, val_range);
X_test  = Dataset_Y(:, :, :, test_range);

% --- 3. 准备并划分标签数据 Y ---
% 合并时延和多普勒索引为 [4600, 3, 2] 的张量
Labels_All = cat(3, Target_Delay_Index, Target_Doppler_Index);

Y_train = Labels_All(train_range, :, :);
Y_val   = Labels_All(val_range, :, :);
Y_test  = Labels_All(test_range, :, :);


% --- 4. 新增：提取并划分 SNR 标签 ---
% 假设你在生成数据时定义了 Target_SNR 变量
% 如果之前没定义，请确保在生成循环中加入了 Target_SNR(i_data) = current_SNR;
SNR_train = Target_SNR(train_range);
SNR_val   = Target_SNR(val_range);
SNR_test  = Target_SNR(test_range);


% --- 5. 打印划分状态确认 ---
fprintf('划分完成（顺序截取）：\n');
fprintf('  - 训练集: 1-4000 (SNR 范围: %d-%d dB)\n', min(SNR_train), max(SNR_train));
fprintf('  - 测试集: 4401-4600 (已提取对应 SNR 标签)\n');

% --- 6. 保存数据 (包含 SNR 标签) ---
save('OTFS_Dataset_Sequential.mat', ...
    'X_train', 'Y_train', 'SNR_train', ...
    'X_val',   'Y_val',   'SNR_val',   ...
    'X_test',  'Y_test',  'SNR_test',  '-v7.3');

fprintf('数据集及 SNR 标签已保存至 OTFS_Dataset_Sequential.mat\n');

%% 打印第一个样本信息
% first_sample_delays   = Target_Delay_Index(1, :);   % 或者 Labels_All(1, :, 1)
% first_sample_dopplers = Target_Doppler_Index(1, :); % 或者 Labels_All(1, :, 2)
% 
% fprintf('\n第一个样本的三个目标参数：\n');
% for i = 1:3
%     fprintf('  目标 %d: [时延索引 = %d, 多普勒索引 = %d]\n', ...
%             i, first_sample_delays(i), first_sample_dopplers(i));
% end
% 
% % 顺便打印对应的实部和虚部最大值（响应您的前一个需求）
% first_sample_data = Dataset_Y(:, :, :, 1);
% max_real = max(reshape(first_sample_data(:, :, 1), [], 1));
% max_imag = max(reshape(first_sample_data(:, :, 2), [], 1));
% 
% fprintf('\n第一个样本的信号强度：\n');
% fprintf('  - 实部最大值: %f\n', max_real);
% fprintf('  - 虚部最大值: %f\n', max_imag);