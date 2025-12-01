% 画图
% --- 1. 初始化环境 ---
clear;      % 清除工作区中的所有变量
clc;        % 清除命令窗口
close all;  % 关闭所有已打开的图形窗口

%% --- 2. 定义文件和参数 ---

SNR=10:2:16;
%%  fig 1
matFiles = {'NMSE1.mat','NMSE2.mat','NMSE3.mat'};
% % 为每个数据集定义图例中显示的标签
legendLabels = {'方案1,QPSK', '方案2,BPSK', '方案3,BPSK+QPSK'};
markerStyles = {'-o', '-s', '-<','--*','--*','--<'};

%%  fig 2
% matFiles = {'2_1622.mat','31_1622.mat','30.75_1622.mat'};
% % % 为每个数据集定义图例中显示的标签
% legendLabels = {'方案2,BPSK', '方案3,BPSK+QPSK','方案3,BPSK+QPSK,0.75'};
% markerStyles = {'-o', '-s', '-<','--*','--*','--<'};


% --- 3. 创建图形并准备绘图 ---
figure;  % 创建一个新的图形窗口

% --- 4. 循环加载数据并绘图 ---
for ifigure = 1:length(matFiles)

    load(matFiles{ifigure});
    fprintf('--- 循环第 %d 次 ---\n', ifigure);
    fprintf('加载文件: %s\n', matFiles{ifigure});
    fprintf('将要使用的图例是: %s\n', legendLabels{ifigure});
    fprintf('将要使用的标记是: %s\n\n', markerStyles{ifigure});
    semilogy(SNR, chanel_err,markerStyles{ifigure} , 'DisplayName', legendLabels{ifigure}, 'LineWidth', 1.5);
    hold on;

    % 清除本次加载的变量，以防对下一次循环造成干扰
    clear ber;
end

% --- 5. 美化和完善图表 ---
hold off; % 绘制完成，结束保持图形状态

% 添加图表标题
title('不同方案下的NMSE');

% 添加 x 轴和 y 轴的标签
xlabel('SNR{(dB)}');
ylabel('NMSE');

% 'Location', 'northeast' 表示图例显示在右上角，您也可以改为 'best' 让MATLAB自动选择最佳位置。
legend('show', 'Location', 'best');

% 显示网格线
grid on;

% % 设置坐标轴的边框样式
% box on;
% xticks(SNR);
% % 2. 定义您想显示的新标签
% new_labels = {'9', '11', '13', '15','17'};
% % % 3. 应用新标签
% xticklabels(new_labels);

