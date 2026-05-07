% 清除工作区和关闭现有图形
clear; clc; close all;

% ==========================================
% 1. 定义参数
% ==========================================
N = 20;
k_vi = 6;          % 整数多普勒频移
kappa_vi = 0.3;    % 分数多普勒频移

% 定义自变量（连续域和离散采样点）
x_cont = linspace(0, 19.5, 2000); 
k_disc = 0:19;                    

% ==========================================
% 2. 计算窗口响应 (调用底部的自定义函数)
% ==========================================
% 整数多普勒
y_int_cont = dirichlet_kernel(x_cont, k_vi, N);
y_int_disc = dirichlet_kernel(k_disc, k_vi, N);

% 分数多普勒
y_frac_cont = dirichlet_kernel(x_cont, k_vi + kappa_vi, N);
y_frac_disc = dirichlet_kernel(k_disc, k_vi + kappa_vi, N);

% ==========================================
% 3. 绘图与美化
% ==========================================
% 设置字体为宋体 (SimSun)
font_name = '宋体'; 
fig = figure('Name', '多普勒响应', 'Position', [100, 100, 700, 500], 'Color', 'w');

% 【核心修改 1】：手动将坐标系撑满，大幅压缩四周留白 [左边距, 下边距, 宽度, 高度]
% 原值为 [0.12 0.12 0.83 0.78]，现修改为极其紧凑的比例
ax = axes('Position', [0.07 0.11 0.91 0.87]); 
hold on; box on;

% 绘制曲线和离散点
p1 = plot(x_cont, y_int_cont, 'r-', 'LineWidth', 1.5);
p2 = plot(k_disc, y_int_disc, 'rx', 'MarkerSize', 6, 'LineWidth', 1.5);
p3 = plot(x_cont, y_frac_cont, 'b--', 'LineWidth', 1.5);
p4 = plot(k_disc, y_frac_disc, 'bo', 'MarkerSize', 6, 'LineWidth', 1.5);

% 设置坐标轴限制和刻度 (顺便将数字设为了更专业的新罗马字体)
xlim([0, 19.5]);
ylim([0, 1.05]);
set(ax, 'XTick', 0:2:18, 'YTick', 0:0.1:1, 'FontSize', 11, 'FontName', 'Times New Roman');

% 设置标签 
xlabel('多普勒频移 \it{k}', 'FontSize', 13, 'FontName', font_name);
ylabel('多普勒域归一化窗响应', 'FontSize', 13, 'FontName', font_name);

% 添加图例 
lgd = legend([p1, p2, p3, p4], ...
    {'整数多普勒(抽头为6)', ...
     '整数 k 处的采样', ...
     '分数多普勒（抽头为6.3）', ...
     '整数 k 处的采样'}, ...
    'FontSize', 11, 'FontName', font_name, 'Location', 'northeast');
lgd.Position = [0.66 0.68 0.25 0.18]; % 微调图例位置以适应更大的画幅

% ==========================================
% 4. 添加精确注释 (分离文本和箭头以避免重叠)
% ==========================================
drawnow; % 确保坐标系已更新以进行准确的坐标转换

% 定义匿名函数：将数据坐标转换为 Figure 归一化坐标
x2n = @(x) ax.Position(1) + ax.Position(3) * (x - ax.XLim(1)) / diff(ax.XLim);
y2n = @(y) ax.Position(2) + ax.Position(4) * (y - ax.YLim(1)) / diff(ax.YLim);

% 4.1 绘制 k=6 处的双向箭头 (功率损失间隙)
y_peak_frac = dirichlet_kernel(6, k_vi + kappa_vi, N);
annotation('doublearrow', [x2n(6) x2n(6)], [y2n(1) y2n(y_peak_frac)], ...
    'Head1Style', 'vback2', 'Head2Style', 'vback2', ...
    'Head1Length', 6, 'Head1Width', 6, 'Head2Length', 6, 'Head2Width', 6);

% 4.2 功率损失文本和指向箭头
text(4.2, 0.88, {'分数多普勒', '引起的', '功率损失'}, 'FontSize', 12, 'FontName', font_name, 'HorizontalAlignment', 'right');
annotation('arrow', [x2n(4.4) x2n(5.8)], [y2n(0.88) y2n(0.92)], ...
    'HeadStyle', 'vback2', 'HeadLength', 6, 'HeadWidth', 6);

% 4.3 功率泄漏文本和两个指向箭头
text(10.5, 0.62, {'分数多普勒引起的功率泄漏'}, 'FontSize', 12, 'FontName', font_name, 'HorizontalAlignment', 'center');
annotation('arrow', [x2n(10.5) x2n(7.2)], [y2n(0.58) y2n(0.38)], ...
    'HeadStyle', 'vback2', 'HeadLength', 6, 'HeadWidth', 6);
annotation('arrow', [x2n(10.5) x2n(9.2)], [y2n(0.58) y2n(0.12)], ...
    'HeadStyle', 'vback2', 'HeadLength', 6, 'HeadWidth', 6);

% ==========================================
% 5. 导出紧凑高清图
% ==========================================
% 【核心修改 2】：导出时自动裁剪绝对边框，生成完美比例图片
% exportgraphics(fig, '多普勒窗响应_无留白版.png', 'Resolution', 600);

% ==========================================
% 本地函数：计算归一化狄利克雷核 (周期性Sinc)
% ==========================================
function y = dirichlet_kernel(x, shift, N)
    arg = x - shift;
    num = sin(pi * arg);
    den = N * sin(pi * arg / N);
    
    y = zeros(size(x));
    tol = 1e-8;
    % 处理分母趋近于0的情况 (极限为1)
    idx_zero = abs(sin(pi * arg / N)) < tol; 
    
    y(~idx_zero) = abs(num(~idx_zero) ./ den(~idx_zero));
    y(idx_zero) = 1; 
end