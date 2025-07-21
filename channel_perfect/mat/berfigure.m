% 画图
% --- 1. 初始化环境 ---
clear;      % 清除工作区中的所有变量
clc;        % 清除命令窗口
close all;  % 关闭所有已打开的图形窗口

%% --- 2. 定义文件和参数 ---

SNR_dB=10:2:16;
%%  fig 1
% matFiles = {'D:\研究生\科研\frac_doppler_channel_est\channel_perfect\mat\ber1qpsk10_16.mat', 'D:\研究生\科研\frac_doppler_channel_est\channel_perfect\mat\ber2bpsk_10_16.mat',...
%  'D:\研究生\科研\frac_doppler_channel_est\channel_perfect\mat\ber3_guard4_cen2\cen2gua4.mat'};
% % 为每个数据集定义图例中显示的标签
% legendLabels = {'方案1,QPSK', '方案2,BPSK', '方案3,BPSK+QPSK'};
%% fig 2
matFiles = {'D:\研究生\科研\frac_doppler_channel_est\channel_perfect\mat\ber3cen2guard4sure.mat', 'D:\研究生\科研\frac_doppler_channel_est\channel_perfect\mat\ber3_guard4_cen2\imori.mat',...
};
% 为每个数据集定义图例中显示的标签
legendLabels = {'方案3,IMsure', '方案3,IM'};
%% fig 3



markerStyles = {'-o', '-s', '-<','--*','--*','--<'}; 
% --- 3. 创建图形并准备绘图 ---
figure;  % 创建一个新的图形窗口


% --- 4. 循环加载数据并绘图 ---
for ifigure = 1:length(matFiles)
    
    % 加载 .mat 文件。文件中的变量将被直接加载到当前工作区。
    % 我们假设每个文件都包含一个名为 'chanel_err' 的变量。
    load(matFiles{ifigure}); 
    
    % 检查变量 'chanel_err' 是否存在于工作区中
    if exist('ber', 'var')
        
        % 检查数据长度是否与 SNR_dB 向量的长度匹配
        if numel(ber) == numel(SNR_dB)
            % 绘制曲线，并使用不同的标记样式以作区分
            % '-o' 表示实线连接圆形标记
            % '-s' 表示实线连接方形标记
            % '-^' 表示实线连接上三角形标记          
            semilogy(SNR_dB, ber,markerStyles{ifigure} , 'DisplayName', legendLabels{ifigure}, 'LineWidth', 1.5);
            hold on;
        else
            % 如果长度不匹配，在命令窗口显示警告信息
            warning('文件 %s 中的 "ber" 数据长度 (%d) 与 SNR_dB 长度 (%d) 不匹配。', ...
                    matFiles{ifigure}, numel(ber), numel(SNR_dB));
        end
        
        % 清除本次加载的变量，以防对下一次循环造成干扰
        clear ber; 
        
    else
        % 如果变量不存在，也在命令窗口显示警告
        warning('在文件 %s 中找不到名为 "ber" 的变量。', matFiles{ifigure});
    end
end

% --- 5. 美化和完善图表 ---
hold off; % 绘制完成，结束保持图形状态

% 添加图表标题
title('不同方案下的BER');

% 添加 x 轴和 y 轴的标签
xlabel('SNR{(dB)}');
ylabel('比特误码率 (BER)');

% 'Location', 'northeast' 表示图例显示在右上角，您也可以改为 'best' 让MATLAB自动选择最佳位置。
legend('show', 'Location', 'best');

% 显示网格线
grid on;

% 设置坐标轴的边框样式
box on;
xticks(SNR_dB); 
% 2. 定义您想显示的新标签
new_labels = {'9', '11', '13', '15'}; 
% 3. 应用新标签
xticklabels(new_labels);

