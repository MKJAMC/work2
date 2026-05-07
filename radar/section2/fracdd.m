% 清除环境
clear; clc; close all;

% 1. 定义 DD 域网格大小
M = 64; 
N = 32; 
[doppler_grid, delay_grid] = meshgrid(1:N, 1:M);

% 2. 模拟三条路径
paths = [
    12.3,  8.2,  0.38;  % 左侧峰
    32.2, 22.3,  0.42;  % 右侧峰 (最高)
    48.7, 18.8,  0.22   % 中间小峰 (最低)
];

% 3. 计算 DD 域响应
Z = zeros(M, N);
for i = 1:size(paths, 1)
    tau = paths(i, 1);
    nu  = paths(i, 2);
    amp = paths(i, 3);
    Z = Z + amp * abs(sinc(delay_grid - tau) .* sinc(doppler_grid - nu));
end

% 4. 绘制 3D 柱状图
figure('Color', 'w', 'Position', [100, 100, 800, 600]); 
h = bar3(Z);

% 5. 调整柱子的颜色和边缘
for k = 1:length(h)
    zdata = h(k).ZData;
    h(k).CData = zdata;             
    h(k).FaceColor = 'interp';      
    h(k).EdgeColor = [0.05 0.05 0.2]; % 边缘线稍微再调暗一点，配合深色底部
    h(k).LineWidth = 0.2;           
end

% 6. 设置颜色映射 (核心修改：自定义非线性细腻蓝色调)
% 定义颜色锚点 (RGB)
colors = [
    0.00, 0.00, 0.25;  % 1. 极暗夜蓝 (最底部的基底)
    0.00, 0.08, 0.50;  % 2. 深海蓝 
    0.00, 0.20, 0.75;  % 3. 中深蓝
    0.00, 0.40, 1.00;  % 4. 亮纯蓝
    0.00, 1.00, 1.00;  % 5. 青色
    0.50, 1.00, 0.40;  % 6. 绿色
    1.00, 1.00, 0.00;  % 7. 黄色
    1.00, 0.10, 0.00;  % 8. 橘红色
    0.60, 0.00, 0.00   % 9. 深红色 (最高峰)
];

% 定义这些颜色锚点在 0~1 范围内的位置 (非线性分布)
% 前面 4 个点距离很近 (0 到 0.25)，强制拉长了底部蓝色的渐变层次！
color_positions = [0, 0.04, 0.12, 0.25, 0.40, 0.60, 0.75, 0.90, 1.00];

% 使用插值生成 256 阶平滑色带
map_custom = interp1(color_positions, colors, linspace(0, 1, 256));
colormap(map_custom);

set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

% B. 单独对标签进行设置，强制指定中文用“宋体”
% 如果 'SimSun' 还是显示方块，请将其替换为 '宋体' (直接写中文名称)
xlabel('多普勒抽头', 'FontName', 'SimSun', 'FontSize', 12, 'FontWeight', 'normal');
ylabel('时延抽头', 'FontName', 'SimSun', 'FontSize', 12, 'FontWeight', 'normal');
zlabel('幅值', 'FontName', 'SimSun', 'FontSize', 12, 'FontWeight', 'normal');

% C. 视角调整
view(-40, 28);
axis tight;
grid on;

