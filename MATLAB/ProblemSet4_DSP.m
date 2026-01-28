% Taha Akhlaq - DSP Problem Set 4
clear; close all; clc;

%% Problem 6

%% Part A

% Haar Filters
h0_haar = [1 1] / sqrt(2);
h1_haar = [1 -1] / sqrt(2);
f0_haar = h0_haar; 
f1_haar = [-1 1] / sqrt(2);

% Check PR condition
t_haar = 0.5 * (conv(f0_haar, h0_haar) + conv(f1_haar, h1_haar));
h0_haar_alt = h0_haar .* [1 -1];
h1_haar_alt = h1_haar .* [1 -1];
a_haar = 0.5 * (conv(f0_haar, h0_haar_alt) + conv(f1_haar, h1_haar_alt));

fprintf('Part a: Haar Wavelet Analysis\n');
fprintf('  H0: [%.3f %.3f]\n', h0_haar);
fprintf('  T(z): '); disp(t_haar);
fprintf('  A(z): '); disp(a_haar);

%% Part B
orders = [4, 8]; 

for i = 1:length(orders)
    N = orders(i);

    % Generate Filters
    [h0, h1, f0, f1] = get_daubechies_filters(N);
    
    h0_flip = fliplr(h0);
    alternating_signs = (-1).^(0:length(h0)-1);
    h1_derived = h0_flip .* alternating_signs;
    
    f0_derived = fliplr(h0);
    f1_derived = fliplr(h1);
    
    err_h1 = min(norm(h1 - h1_derived), norm(h1 + h1_derived));
    err_f0 = min(norm(f0 - f0_derived), norm(f0 + f0_derived));
    
    fprintf('\n');
    fprintf('Part b: Analysis for N = %d\n', N);
    fprintf('1. Filter Relation Errors:\n');
    fprintf('   H1 error: %e\n', err_h1);
    fprintf('   F0 error: %e\n', err_f0);
    
    % Plots
    num_points = 10000;
    [H0_w, w] = freqz(h0, 1, num_points, 'whole');
    [H1_w, ~] = freqz(h1, 1, num_points, 'whole');
    
    w_half = w(1:num_points/2);
    H0_mag = abs(H0_w(1:num_points/2));
    H1_mag = abs(H1_w(1:num_points/2));
    
    figure;
    plot(w_half, H0_mag, 'b', 'LineWidth', 1.5); hold on;
    plot(w_half, H1_mag, 'r', 'LineWidth', 1.5);
    grid on;
    title(['Magnitude Responses N=', num2str(N)]);
    xlabel('\omega (rad)'); ylabel('Magnitude');
    legend('|H_0(\omega)|', '|H_1(\omega)|');
    xlim([0 pi]);
    
    % Power Check
    power_sum = H0_mag.^2 + H1_mag.^2;
    max_power_err = max(abs(power_sum - 2));
    fprintf('3. Power Error: %e\n', max_power_err);
    
    % Polyphase
    e00 = h0(1:2:end);
    e01 = h0(2:2:end);
    e10 = h1(1:2:end);
    e11 = h1(2:2:end);
    
    r00 = fliplr(e00); 
    r01 = fliplr(e10); 
    r10 = fliplr(e01); 
    r11 = fliplr(e11); 
    
    % Compute P(z) = R(z)E(z)
    p00 = conv(r00, e00) + conv(r01, e10);
    p01 = conv(r00, e01) + conv(r01, e11);
    p10 = conv(r10, e00) + conv(r11, e10);
    p11 = conv(r10, e01) + conv(r11, e11);
    
    % Check Errors
    off_diag_err = max([max(abs(p01)), max(abs(p10))]);
    [max_p00, idx_p00] = max(abs(p00));
    diag_err = sum(abs(p00)) - max_p00; 
    
    d = p00(idx_p00);
    M_delay = idx_p00 - 1; 
    
    fprintf('6. P(z) = R(z)E(z):\n');
    fprintf('   Constant: %.4f\n', d);
    fprintf('   Delay: %d\n', M_delay);
    fprintf('   Max Error: %e\n', max(off_diag_err, diag_err));
    
    % PR Check
    alt_seq_h0 = (-1).^(0:length(h0)-1);
    alt_seq_h1 = (-1).^(0:length(h1)-1);
    h0_minus = h0 .* alt_seq_h0;
    h1_minus = h1 .* alt_seq_h1;
    
    T_total = 0.5 * (conv(f0, h0) + conv(f1, h1));
    A_total = 0.5 * (conv(f0, h0_minus) + conv(f1, h1_minus));
    
    [max_T, idx_T] = max(abs(T_total));
    T_err = sum(abs(T_total)) - max_T;
    A_err = max(abs(A_total));
    
    fprintf('7. PR Check:\n');
    fprintf('   T(z) Error: %e\n', T_err);
    fprintf('   A(z) Error: %e\n', A_err);
    
    H0_sq = H0_mag.^2;
    dw = w_half(2) - w_half(1);
    d1 = diff(H0_sq) / dw; d1_plot = [d1; 0];
    d2 = diff(d1) / dw;    d2_plot = [d2; 0; 0];
    
    figure;
    subplot(2,1,1); plot(w_half, d1_plot, 'b'); grid on;
    title(['1st Derivative of |H_0|^2 (N=', num2str(N), ')']); xlim([0 pi]);
    subplot(2,1,2); plot(w_half, d2_plot, 'r'); grid on;
    title(['2nd Derivative of |H_0|^2 (N=', num2str(N), ')']); xlim([0 pi]);
end

% Tree Structure

% N=4
[h0, h1, ~, ~] = get_daubechies_filters(4);

% Upsample
upsample_filter = @(h, M) reshape([h; zeros(M-1, length(h))], 1, []);

% Level 2
h0_z2 = upsample_filter(h0, 2);
h1_z2 = upsample_filter(h1, 2);

G0 = conv(h0, h0_z2); 
G1 = conv(h0, h1_z2); 
G2 = conv(h1, h0_z2); 
G3 = conv(h1, h1_z2);

figure;
subplot(2,1,1); hold on;
plot_freq_response(G0, 'b'); plot_freq_response(G1, 'g');
plot_freq_response(G2, 'r'); plot_freq_response(G3, 'k');
title('Magnitude Responses: 2-Level Tree');
xlabel('\omega'); ylabel('|G(\omega)|'); grid on; xlim([0 pi]);

% Level 3
h0_z4 = upsample_filter(h0, 4);
h1_z4 = upsample_filter(h1, 4);

% Compute 8 branches
filters_L3 = cell(1,8);
count = 1;
G_level2 = {G0, G1, G2, G3};
for k = 1:4
    filters_L3{count} = conv(G_level2{k}, h0_z4);
    filters_L3{count+1} = conv(G_level2{k}, h1_z4);
    count = count + 2;
end

subplot(2,1,2); hold on;
colors = lines(8);
for k = 1:8
    plot_freq_response(filters_L3{k}, colors(k,:));
end
title('Magnitude Responses: 3-Level Tree');
xlabel('\omega'); ylabel('|G(\omega)|'); grid on; xlim([0 pi]);


% N=8
[h0, h1, ~, ~] = get_daubechies_filters(8);

% Upsample
upsample_filter = @(h, M) reshape([h; zeros(M-1, length(h))], 1, []);

% 2 Level
h0_z2 = upsample_filter(h0, 2);
h1_z2 = upsample_filter(h1, 2);

G0 = conv(h0, h0_z2); 
G1 = conv(h0, h1_z2); 
G2 = conv(h1, h0_z2); 
G3 = conv(h1, h1_z2);

figure;
subplot(2,1,1); hold on;
plot_freq_response(G0, 'b'); plot_freq_response(G1, 'g');
plot_freq_response(G2, 'r'); plot_freq_response(G3, 'k');
title('Magnitude Responses: 2-Level Tree (N=8)');
xlabel('\omega'); ylabel('|G(\omega)|'); grid on; xlim([0 pi]);

% 3 Level
h0_z4 = upsample_filter(h0, 4);
h1_z4 = upsample_filter(h1, 4);

% Compute 8 branches
filters_L3 = cell(1,8);
count = 1;
G_level2 = {G0, G1, G2, G3};

for k = 1:4
    filters_L3{count} = conv(G_level2{k}, h0_z4);
    filters_L3{count+1} = conv(G_level2{k}, h1_z4);
    count = count + 2;
end

subplot(2,1,2); hold on;
colors = lines(8);
for k = 1:8
    plot_freq_response(filters_L3{k}, colors(k,:));
end
title('Magnitude Responses: 3-Level Tree (N=8)');
xlabel('\omega'); ylabel('|G(\omega)|'); grid on; xlim([0 pi]);

%% Functions

function plot_freq_response(h, color_spec)
    [H, w] = freqz(h, 1, 5000, 'whole');
    w_half = w(1:end/2);
    mag = abs(H(1:end/2));
    if ischar(color_spec)
        plot(w_half, mag, color_spec);
    else
        plot(w_half, mag, 'Color', color_spec);
    end
end

% function to replace wfilters
function [h0, h1, f0, f1] = get_daubechies_filters(N)
    if N == 4
        h0 = [0.230377813309, 0.714846570553, 0.630880767930, -0.027983769417, ...
              -0.187034811719, 0.030841381836, 0.032883011667, -0.010597401785];
    elseif N == 8
        h0 = [0.054415842243, 0.312871590914, 0.675630736297, 0.585354683654, ...
              -0.015829105256, -0.284015542962, 0.000472484574, 0.128747426620, ...
              -0.017369301002, -0.044088253931, 0.013981027917, 0.008746094047, ...
              -0.004870352993, -0.000391740373, 0.000675449406, -0.000117476784];
    else
        error('Invalid N');
    end
    
    L = length(h0);
    h0_flip = fliplr(h0);
    alt_signs = (-1).^(0:L-1);
    h1 = h0_flip .* alt_signs;
    
    f0 = fliplr(h0);
    f1 = fliplr(h1);
end