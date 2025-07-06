% 画图
% --- 1. 初始化环境 ---
clear;      % 清除工作区中的所有变量
clc;        % 清除命令窗口
close all;  % 关闭所有已打开的图形窗口

%% --- 2. 定义文件和参数 ---
SNR_dB = 0:5:15;

% 将您的三个 .mat 文件名存储在一个元胞数组 (cell array) 中
matFiles = {'first_snr_15.mat', 'second_snr_15.mat', 'third_1.mat','third_0.75.mat','third_0.5.mat'};

% 为每个数据集定义图例中显示的标签
legendLabels = {'方案1', '方案2', '方案3,{1}','方案三,{0.75}','方案三,{0.5}'};
markerStyles = {'-o', '-s', '-^','-square','-*'}; 
% --- 3. 创建图形并准备绘图 ---
figure;  % 创建一个新的图形窗口
hold on; % 允许在同一张图上绘制多条曲线

% --- 4. 循环加载数据并绘图 ---
for ifigure = 1:length(matFiles)
    
    % 加载 .mat 文件。文件中的变量将被直接加载到当前工作区。
    % 我们假设每个文件都包含一个名为 'chanel_err' 的变量。
    load(matFiles{ifigure}); 
    
    % 检查变量 'chanel_err' 是否存在于工作区中
    if exist('chanel_err', 'var')
        
        % 检查数据长度是否与 SNR_dB 向量的长度匹配
        if numel(chanel_err) == numel(SNR_dB)
            % 绘制曲线，并使用不同的标记样式以作区分
            % '-o' 表示实线连接圆形标记
            % '-s' 表示实线连接方形标记
            % '-^' 表示实线连接上三角形标记
    
      
            plot(SNR_dB, chanel_err, markerStyles{ifigure}, 'DisplayName', legendLabels{ifigure}, 'LineWidth', 1.5);
            
        else
            % 如果长度不匹配，在命令窗口显示警告信息
            warning('文件 %s 中的 "chanel_err" 数据长度 (%d) 与 SNR_dB 长度 (%d) 不匹配。', ...
                    matFiles{ifigure}, numel(chanel_err), numel(SNR_dB));
        end
        
        % 清除本次加载的变量，以防对下一次循环造成干扰
        clear chanel_err; 
        
    else
        % 如果变量不存在，也在命令窗口显示警告
        warning('在文件 %s 中找不到名为 "chanel_err" 的变量。', matFiles{ifigure});
    end
end

% --- 5. 美化和完善图表 ---
hold off; % 绘制完成，结束保持图形状态

% 添加图表标题
title('不同方案下的NMSE');

% 添加 x 轴和 y 轴的标签
xlabel('SNR_{d}(dB)');
ylabel('NMSE ');
% 显示图例。'DisplayName' 已在 plot 函数中设置好了。
% 'Location', 'northeast' 表示图例显示在右上角，您也可以改为 'best' 让MATLAB自动选择最佳位置。
legend('show', 'Location', 'best');

% 显示网格线
grid on;

% 设置坐标轴的边框样式
box on;

% --- 6. 保存图像 (可选) ---
% 如果需要将生成的图像保存为文件（如PNG格式），可以取消下面这行代码的注释
% saveas(gcf, 'channel_error_plot.png'); % gcf 获取当前图形句柄