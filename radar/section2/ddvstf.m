% 终极对比版：高移动性 vs 零多普勒 (引入 tiledlayout 去除留白，完美排版)
clear; clc; close all;
rng(123); 

%% 1. 系统全局参数设置
M = 64;             
N = 64;             
df = 15e3;          
T = 1/df;           

% TF域保留物理单位
f_vec = linspace(-0.5*M*df, 0.5*M*df, M);
t_vec = linspace(0, N*T, N);
[F, Time] = meshgrid(f_vec, t_vec);

%% ================= 场景一：多普勒和时延均不为 0 (时变信道) =================
% 将所有负数多普勒频移改为正数，确保多普勒抽头 k > 0
taps_tv = [
    1.0,  0.01*T,   0.02/T,  rand()*2*pi;   % 主径
    0.8,  0.04*T,   0.03/T,  rand()*2*pi;   % 反射径1 
    0.6,  0.07*T,   0.05/T,  rand()*2*pi;   % 反射径2
    0.5,  0.09*T,   0.01/T,  rand()*2*pi;   % 反射径3 
    0.3,  0.12*T,   0.06/T,  rand()*2*pi    % 反射径4
];
gains_complex_tv = taps_tv(:, 1) .* exp(1i * taps_tv(:, 4)); 
taus_tv = taps_tv(:, 2); 
nus_tv = taps_tv(:, 3);

% DD域转换为离散抽头 (Tap)
taus_tap_tv = taus_tv * (M / T);  
nus_tap_tv  = nus_tv * (N * T);   

H_ft_tv = zeros(N, M);
for l = 1:size(taps_tv, 1)
    H_ft_tv = H_ft_tv + gains_complex_tv(l) * exp(-1i * 2 * pi * F * taus_tv(l)) .* exp(1i * 2 * pi * nus_tv(l) * Time);
end
H_mag_tv = abs(H_ft_tv);

% --- 绘制第一幅图 ---
figure('Name', '高移动性信道', 'Position', [100, 500, 950, 420], 'Color', 'w');
% 【核心修改 1】：使用 tiledlayout 彻底消除留白
t1 = tiledlayout(1, 2, 'Padding', 'tight', 'TileSpacing', 'tight');

% 图 1-(a)
nexttile; % 替代 subplot
stem3(taus_tap_tv, nus_tap_tv, abs(gains_complex_tv), 'filled', 'MarkerSize', 6, 'LineWidth', 2.5, 'Color', [0.1 0.6 0.3]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
xlabel('时延抽头', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('多普勒抽头', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold'); 
zlabel('幅度', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold');
% 【核心修改 2】：使用 title 并在底部居中显示，替代 annotation
title('(a) 时延-多普勒域 h(l, k)', 'FontName', '宋体', 'FontSize', 14, 'FontWeight', 'normal', 'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
grid on; 
xlim([-1, 9]); 
ylim([-1, 5]); 
zlim([0, 1.2]); 
view(-35, 30);

% 图 1-(b)
nexttile; % 替代 subplot
surf(f_vec, t_vec, H_mag_tv);
shading interp; colormap(jet); 
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
xlabel('频率 f (Hz)', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('时间 t (s)', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold'); 
zlabel('幅度', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold');
title('(b) 时频域 H(f, t)', 'FontName', '宋体', 'FontSize', 14, 'FontWeight', 'normal', 'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
grid on; axis tight; zlim([0, max(H_mag_tv(:))*1.1]); 
view(-35, 30);


%% ================= 场景二：多普勒为 0，时延不为 0 (时不变信道) =================
taps_ti = taps_tv;       
taps_ti(:, 3) = 0;       
gains_complex_ti = taps_ti(:, 1) .* exp(1i * taps_ti(:, 4)); 
taus_ti = taps_ti(:, 2); 
nus_ti = taps_ti(:, 3);

taus_tap_ti = taus_ti * (M / T);  
nus_tap_ti  = nus_ti * (N * T);   

H_ft_ti = zeros(N, M);
for l = 1:size(taps_ti, 1)
    H_ft_ti = H_ft_ti + gains_complex_ti(l) * exp(-1i * 2 * pi * F * taus_ti(l)) .* exp(1i * 2 * pi * nus_ti(l) * Time);
end
H_mag_ti = abs(H_ft_ti);

% --- 绘制第二幅图 ---
figure('Name', '零多普勒信道', 'Position', [150, 450, 950, 420], 'Color', 'w');
% 【核心修改 1】：使用 tiledlayout
t2 = tiledlayout(1, 2, 'Padding', 'tight', 'TileSpacing', 'tight');

% 图 2-(a)
nexttile; 
stem3(taus_tap_ti, nus_tap_ti, abs(gains_complex_ti), 'filled', 'MarkerSize', 6, 'LineWidth', 2.5, 'Color', [0.8 0.1 0.1]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
xlabel('时延抽头', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('多普勒抽头', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold'); 
zlabel('幅度', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold');
title('(a) 时延-多普勒域 h(l, k)', 'FontName', '宋体', 'FontSize', 14, 'FontWeight', 'normal', 'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
grid on; 
xlim([-1, 9]); 
ylim([-1, 5]); 
zlim([0, 1.2]); 
view(-35, 30);

% 图 2-(b)
nexttile; 
surf(f_vec, t_vec, H_mag_ti);
shading interp; colormap(jet);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
xlabel('频率 f (Hz)', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('时间 t (s)', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold'); 
zlabel('幅度', 'FontName', '宋体', 'FontSize', 12, 'FontWeight', 'bold');
title('(b) 时频域 H(f, t)', 'FontName', '宋体', 'FontSize', 14, 'FontWeight', 'normal', 'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
grid on; axis tight; zlim([0, max(H_mag_ti(:))*1.1]); 
view(-35, 30);