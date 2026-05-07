clc; clear all; close all;

% 1. 加载划分后的数据集
load('OTFS_Dataset_Sequential.mat');

% 设定要打印的样本数量
num_samples_to_print = 3;

fprintf('==========================================\n');
fprintf('  OTFS 测试集 - 前 %d 个样本详细信息\n', num_samples_to_print);
fprintf('==========================================\n');

% 2. 循环遍历前 3 个样本
for sample_idx = 1:num_samples_to_print
    
    % --- 数据提取 ---
    
    % 获取该样本对应的原始 SNR
    sample_snr = SNR_test(sample_idx);
    
    % 获取该样本的三个目标标签 [3x2]
    % 第一列是 Delay Index, 第二列是 Doppler Index
    sample_labels = squeeze(Y_test(sample_idx, :, :));
    
    % 获取该样本的输入数据 (M x N x 2)
    sample_data = X_test(:, :, :, sample_idx);
    
    % 计算最大值用于归一化检查
    max_real = max(reshape(sample_data(:,:,1), [], 1));
    max_imag = max(reshape(sample_data(:,:,2), [], 1));
    
    % --- 打印结果 ---
    fprintf('\n------------------------------------------\n');
    fprintf('  样本索引 (Index): %d\n', sample_idx);
    fprintf('------------------------------------------\n');
    fprintf('当前样本信噪比 (SNR): %d dB\n', sample_snr);
    
    fprintf('目标参数 (Ground Truth):\n');
    for i = 1:3
        fprintf('  目标 %d: [时延索引 = %.4f, 多普勒索引 = %.4f]\n', ...
                i, sample_labels(i, 1), sample_labels(i, 2));
    end
    % 
    % fprintf('输入信号强度 (Normalization Check):\n');
    % fprintf('  - 实部 (Real) 最大值: %.6f\n', max_real);
    % fprintf('  - 虚部 (Imag) 最大值: %.6f\n', max_imag);
end

fprintf('\n==========================================\n');
fprintf('  打印结束\n');
fprintf('==========================================\n');