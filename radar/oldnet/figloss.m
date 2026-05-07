clc;clear all;close all;
% 1. 读取 CSV 数据
data = readtable('loss_log.csv');
epochs = data.epoch;
train_loss = data.train_loss;
val_loss = data.val_loss;

% 2. 创建画布
figure('Color', 'w', 'Position', [100, 100, 800, 500]);

% 3. 使用 semilogy 绘制对数纵坐标曲线
% semilogy 会自动将 y 轴设为对数刻度，并处理网格
p1 = semilogy(epochs, train_loss, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0, 0.447, 0.741]); 
hold on;
p2 = semilogy(epochs, val_loss, '-s', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.85, 0.325, 0.098]);

% 4. 图表修饰
grid on;            % 开启主网格
grid minor;         % 开启次网格（对数坐标下非常重要，能辅助看清数值）

xlabel('Epoch', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Loss (Log Scale)', 'FontSize', 12, 'FontWeight', 'bold');
title('Training and Validation Loss Trend', 'FontSize', 14);

% 设置横坐标刻度为整数
set(gca, 'XTick', epochs);

% 设置图例
legend([p1, p2], {'Training Loss', 'Validation Loss'}, 'Location', 'northeast', 'FontSize', 10);

% 5. 自动调整 y 轴范围，使其更美观
ylim([min([train_loss; val_loss])*0.8, max([train_loss; val_loss])*1.2]);

hold off;

% 提示：如果需要保存图片，取消下面一行的注释
% exportgraphics(gcf, 'loss_trend_matlab.png', 'Resolution', 300);