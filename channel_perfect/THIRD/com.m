% --- MATLAB Code to Reproduce Computational Complexity Plot ---

% --- 1. Initialization ---
clear;          % Clear workspace variables
clc;            % Clear command window
close all;      % Close all figures

% --- 2. Define Parameters ---
N = 32;                         % Given parameter from the plot
M = 25:200;     % Define the range for M (200 points for a smooth curve)

% --- 3. Calculate Computational Complexity ---

% Conventional MMSE equalizer complexity is dominated by the inversion of an
% MN x MN matrix, which is on the order of O((MN)^3).
comp_conv = (M * N).^2;

% The proposed low-complexity MMSE equalizer uses FFTs to diagonalize the
% channel matrix. The complexity is dominated by the 2D-FFT, which is on
% the order of O(MN * log(MN)).
% Note: We use log2 as is standard for FFT complexity analysis.
comp_prop = (M * N) .* log2(M * N);

% --- 4. Plotting ---
figure; % Create a new figure window

% Plot using a logarithmic scale for the y-axis (semilogy)
% Plot the conventional MMSE curve in red
semilogy(M, comp_conv, 'r-', 'LineWidth', 2);
hold on; % Hold the plot to add the second curve

% Plot the proposed low-complexity MMSE curve in blue
semilogy(M, comp_prop, 'b-', 'LineWidth', 2);

% --- 5. Formatting the Plot ---
hold off; % Release the plot
grid on;  % Turn on the grid

% Set the axis limits to match the source image
axis([25 200 1e4 1e9]); 

% Add labels and title
% xlabel('M', 'FontSize', 12);
% ylabel('No. operations', 'FontSize', 12);
% title('Computational complexity of conventional and proposed MMSE equalizers', 'FontSize', 12);

ylabel('运算次数', 'FontSize', 12);
title('计算复杂度对比', 'FontSize', 12);
% Add the legend
% legend('Conv. MMSE equal.', 'Prop. low comp. MMSE equal.', ...
%        'Location', 'northwest', 'FontSize', 10);
legend('传统MMSE', '低复杂度MMSE', ...
       'Location', 'northwest', 'FontSize', 10);

% Add the text annotation inside the plot
text(150, 10^7, 'N = 32', 'FontSize', 12, 'HorizontalAlignment', 'center');

% Improve the overall appearance
set(gca, 'FontSize', 10); % Set the font size for the axes tick labels