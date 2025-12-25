clc;clear all;close all;
% 1. 加载划分后的数据集
load('OTFS_Dataset_Sequential.mat');

% 2. 提取测试集的第一个样本 (Index 1 in X_test/Y_test)
% X_test 维度: [M, N, 2, 200]
% Y_test 维度: [200, 3, 2]
% SNR_test 维度: [200, 1]

sample_idx = 1; % 测试集中的第 1 个

% 获取该样本对应的原始 SNR
sample_snr = SNR_test(sample_idx);

% 获取该样本的三个目标标签 [3x2]
% 第一列是 Delay Index, 第二列是 Doppler Index
sample_labels = squeeze(Y_test(sample_idx, :, :));

% 获取该样本的输入数据 (M x N x 2)
sample_data = X_test(:, :, :, sample_idx);
max_real = max(reshape(sample_data(:,:,1), [], 1));
max_imag = max(reshape(sample_data(:,:,2), [], 1));

% --- 打印结果 ---
fprintf('==========================================\n');
fprintf('  OTFS 测试集 - 第 1 个样本详细信息\n');
fprintf('==========================================\n');
fprintf('当前样本信噪比 (SNR): %d dB\n', sample_snr);
fprintf('\n目标参数 (Ground Truth):\n');
for i = 1:3
    fprintf('  目标 %d: [时延索引 = %.4f, 多普勒索引 = %.4f]\n', ...
            i, sample_labels(i, 1), sample_labels(i, 2));
end

fprintf('\n输入信号强度 (Normalization Check):\n');
fprintf('  - 实部 (Real) 最大值: %.6f\n', max_real);
fprintf('  - 虚部 (Imag) 最大值: %.6f\n', max_imag);
fprintf('==========================================\n');