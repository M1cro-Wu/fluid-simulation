clear; clc; close all;
%% load the data
load('Data_Fig1.mat'); 

%% Identify Crests, Troughs, and Zero-Crossings
[crests, t_crests] = findpeaks(eta_max_wave, t_max_wave);
[troughs_inv, t_troughs] = findpeaks(-eta_max_wave, t_max_wave);
troughs = -troughs_inv;

idx_zero = find(eta_max_wave(1:end-1) .* eta_max_wave(2:end) < 0);
t_zero = t_max_wave(idx_zero) - eta_max_wave(idx_zero) .* (t_max_wave(idx_zero+1) - t_max_wave(idx_zero)) ./ (eta_max_wave(idx_zero+1) - eta_max_wave(idx_zero));

is_up_crossing = eta_max_wave(idx_zero) < 0; 
t_z_up = t_zero(is_up_crossing);
t_z_down = t_zero(~is_up_crossing);

%% Extract Specific Points Around the Main Wave
[max_crest, idx_max] = max(crests);
t_c = t_crests(idx_max);

% Adjacent Crests
t_c_prev = t_crests(idx_max - 1); c_prev = crests(idx_max - 1);
t_c_next = t_crests(idx_max + 1); c_next = crests(idx_max + 1);

% Troughs immediately bounding the main crest
t_t_prev = max(t_troughs(t_troughs < t_c)); t_prev = troughs(t_troughs == t_t_prev);
t_t_next = min(t_troughs(t_troughs > t_c)); t_next = troughs(t_troughs == t_t_next);

% Zero crossings that strictly encompass the main crest
t_z_up_start = max(t_z_up(t_z_up < t_c)); 
t_z_up_end = min(t_z_up(t_z_up > t_c)); 

t_z_down_start = max(t_z_down(t_z_down < t_c)); 
t_z_down_end = min(t_z_down(t_z_down > t_c));   

%% Evaluate T and H
% Evaluate Periods
T1 = t_t_next - t_t_prev;                 % Trough to Trough
T2 = t_z_down_end - t_z_down_start;       % Zero-Down to Zero-Down
T3 = t_z_up_end - t_z_up_start;           % Zero-Up to Zero-Up
T4 = t_c - t_c_prev;                      % Pre-crest to Main Crest
T5 = t_c_next - t_c;                      % Main Crest to After Crest

% Evaluate Heights (Max - Min within their respective T interval)
H1 = max_crest - min(t_prev, t_next);     % Max crest - lowest bounding trough
H2 = max_crest - t_prev;                  % Max crest - preceding trough
H3 = max_crest - t_next;                  % Max crest - succeeding trough
H4 = max_crest - t_prev;                  % Max crest - preceding trough
H5 = max_crest - t_next;                  % Max crest - succeeding trough

% Store in Arrays
T = [T1, T2, T3, T4, T5];
H = [H1, H2, H3, H4, H5];

%% Plot for H and T
figure('Position', [100, 100, 1000, 600]);
plot(t_max_wave, eta_max_wave, 'b-', 'LineWidth', 1.2); hold on;
plot(t_crests, crests, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r'); 
plot(t_troughs, troughs, 'yo', 'MarkerSize', 5, 'MarkerFaceColor', 'y');
plot(t_zero, zeros(size(t_zero)), 'mo', 'MarkerSize', 5, 'MarkerFaceColor', 'm');

% ==========================================
% --- Draw T Dimensions (Periods) ---
% ==========================================

% T1: Trough to Trough (includes main crest)
plot([t_t_prev, t_t_next], [-3.8 -3.8], 'k-', 'LineWidth', 1.5);
plot([t_t_prev t_t_prev], [t_prev -4.2], 'k--', 'LineWidth', 1);
plot([t_t_next t_t_next], [t_next -4.2], 'k--', 'LineWidth', 1);
text(mean([t_t_prev, t_t_next]), -3.5, 'T_1', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

% T2: Zero-Down to Zero-Down (includes main crest & preceding trough)
plot([t_z_down_start, t_z_down_end], [-0.8 -0.8], 'k-', 'LineWidth', 1.5);
plot([t_z_down_start t_z_down_start], [0 -1.1], 'k--', 'LineWidth', 1);
plot([t_z_down_end t_z_down_end], [0 -1.1], 'k--', 'LineWidth', 1);
text(mean([t_z_down_start, t_z_down_end]), -1.2, 'T_2', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

% T3: Zero-Up to Zero-Up (includes main crest & succeeding trough)
plot([t_z_up_start, t_z_up_end], [0.8 0.8], 'k-', 'LineWidth', 1.5);
plot([t_z_up_start t_z_up_start], [0 1.1], 'k--', 'LineWidth', 1);
plot([t_z_up_end t_z_up_end], [0 1.1], 'k--', 'LineWidth', 1);
text(mean([t_z_up_start, t_z_up_end]), 1.2, 'T_3', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

% T4: Pre-crest to Main Crest
plot([t_c_prev, t_c], [6.8 6.8], 'k--', 'LineWidth', 1.2);
plot([t_c_prev t_c_prev], [c_prev 7.2], 'k--', 'LineWidth', 1);
plot([t_c t_c], [max_crest 7.2], 'k--', 'LineWidth', 1);
text(mean([t_c_prev, t_c]), 7.2, 'T_4', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

% T5: Main Crest to After Crest
plot([t_c, t_c_next], [6.8 6.8], 'k--', 'LineWidth', 1.2);
plot([t_c_next t_c_next], [c_next 7.2], 'k--', 'LineWidth', 1);
text(mean([t_c, t_c_next]), 7.2, 'T_5', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');


% ==========================================
% --- Draw H Dimensions (Heights) ---
% All heights defined as Max - Min within their respective T interval
% ==========================================

% H1: Max - Min within T1 (Trough to Trough)
min_h1 = min(t_prev, t_next);
h1_x = t_c + 1.0; % Offset slightly right
plot([h1_x, h1_x], [min_h1, max_crest], 'k-', 'LineWidth', 1.2);
plot([h1_x-0.5 h1_x+0.5], [max_crest max_crest], 'k-', 'LineWidth', 1); % Top cap
plot([h1_x-0.5 h1_x+0.5], [min_h1 min_h1], 'k-', 'LineWidth', 1);       % Bottom cap
text(h1_x+0.3, mean([min_h1, max_crest]), 'H_1', 'HorizontalAlignment', 'left', 'FontSize', 12, 'FontWeight', 'bold');

% H2: Max - Min within T2 (Zero-Down to Zero-Down) -> Min is preceding trough
min_h2 = t_prev;
h2_x = t_c - 1.0; % Offset slightly left
plot([h2_x, h2_x], [min_h2, max_crest], 'k-', 'LineWidth', 1.2);
plot([h2_x-0.5 h2_x+0.5], [max_crest max_crest], 'k-', 'LineWidth', 1); % Top cap
plot([h2_x-0.5 h2_x+0.5], [min_h2 min_h2], 'k-', 'LineWidth', 1);       % Bottom cap
text(h2_x-0.3, mean([min_h2, max_crest]), 'H_2', 'HorizontalAlignment', 'right', 'FontSize', 12, 'FontWeight', 'bold');

% H3: Max - Min within T3 (Zero-Up to Zero-Up) -> Min is succeeding trough
min_h3 = t_next;
h3_x = t_c + 2.5; % Offset further right
plot([h3_x, h3_x], [min_h3, max_crest], 'k-', 'LineWidth', 1.2);
plot([h3_x-0.5 h3_x+0.5], [max_crest max_crest], 'k-', 'LineWidth', 1); % Top cap
plot([h3_x-0.5 h3_x+0.5], [min_h3 min_h3], 'k-', 'LineWidth', 1);       % Bottom cap
text(h3_x+0.3, mean([min_h3, max_crest]), 'H_3', 'HorizontalAlignment', 'left', 'FontSize', 12, 'FontWeight', 'bold');

% H4: Max - Min within T4 (Pre-crest to Main crest) -> Min is preceding trough
min_h4 = t_prev;
h4_x = t_c - 2.5; % Offset further left
plot([h4_x, h4_x], [min_h4, max_crest], 'k-', 'LineWidth', 1.2);
plot([h4_x-0.5 h4_x+0.5], [max_crest max_crest], 'k-', 'LineWidth', 1); % Top cap
plot([h4_x-0.5 h4_x+0.5], [min_h4 min_h4], 'k-', 'LineWidth', 1);       % Bottom cap
text(h4_x-0.3, mean([min_h4, max_crest]), 'H_4', 'HorizontalAlignment', 'right', 'FontSize', 12, 'FontWeight', 'bold');

% H5: Max - Min within T5 (Main crest to After crest) -> Min is succeeding trough
min_h5 = t_next;
h5_x = t_c + 4.0; % Offset furthest right
plot([h5_x, h5_x], [min_h5, max_crest], 'k-', 'LineWidth', 1.2);
plot([h5_x-0.5 h5_x+0.5], [max_crest max_crest], 'k-', 'LineWidth', 1); % Top cap
plot([h5_x-0.5 h5_x+0.5], [min_h5 min_h5], 'k-', 'LineWidth', 1);       % Bottom cap
text(h5_x+0.3, mean([min_h5, max_crest]), 'H_5', 'HorizontalAlignment', 'left', 'FontSize', 12, 'FontWeight', 'bold');

% Graph formatting
xlim([-12 10]); % Zoom in to clearly see all labels around the crest
ylim([-5 8.5]); 
grid on;
xlabel('Time [s]');
ylabel('Eta Surface Elevation [m]');
title('Definition of Wave Periods and Heights around the Largest Crest');
hold off;

