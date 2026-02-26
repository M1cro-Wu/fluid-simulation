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
%% validity diagram

% =========================================================
% Le Méhauté (1976) Validity Diagram Overlay (With Legend)
% =========================================================

d = 30;           % Water depth [m]
g = 9.81;         % Acceleration due to gravity [m/s^2]

% Dimensionless parameters for the evaluated waves
x_dim = d ./ (g * T.^2);
y_dim = H ./ (g * T.^2);

figure('Position', [150, 150, 900, 700]);
hold on;

% =========================================================
% 1. Construct the Le Méhauté Background Boundaries
% =========================================================
% Create a dense array for the X-axis to draw smooth curves
x_bg = logspace(-3, 0, 500);

% --- Breaking Limit Curve ---
y_break = min(0.78 .* x_bg, 0.027);
plot(x_bg, y_break, 'k-', 'LineWidth', 2.5, 'HandleVisibility', 'off');

% --- Shallow & Deep Water Vertical Boundaries ---
xline(0.0025, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0.08, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

% --- Ursell Parameter Boundary (Stokes vs. Cnoidal) ---
y_ursell = 26 .* (x_bg.^2);
y_ursell(y_ursell > y_break) = NaN; 
plot(x_bg, y_ursell, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% --- Stokes Order Empirical Boundaries ---
y_stokes_high = 0.1 .* (x_bg.^1.5); 
y_stokes_high(y_stokes_high > y_break) = NaN;
plot(x_bg, y_stokes_high, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');

y_stokes_low = 0.005 .* (x_bg.^1.2); 
y_stokes_low(y_stokes_low > y_break) = NaN;
plot(x_bg, y_stokes_low, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% --- Add Text Labels for the Regions ---
text(0.0012, 0.04, 'BREAKING WAVES', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
text(0.0012, 0.001, 'CNOIDAL / STREAM', 'FontSize', 10, 'FontWeight', 'bold');
text(0.01, 0.005, 'STOKES 3rd / 5th', 'FontSize', 10, 'FontWeight', 'bold');
text(0.02, 0.0003, 'STOKES 2nd', 'FontSize', 10, 'FontWeight', 'bold');
text(0.15, 0.00003, 'LINEAR THEORY', 'FontSize', 10, 'FontWeight', 'bold');

text(0.0011, 1.5e-5, 'SHALLOW', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
text(0.01, 1.5e-5, 'TRANSITIONAL', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
text(0.2, 1.5e-5, 'DEEP', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);


% =========================================================
% 2. Overlay Evaluated T and H Points with Projection Lines
% =========================================================
colors = {'b', 'r', [0 0.5 0], 'm', [0.85 0.33 0.1]};

for i = 1:5
    % Create label string for the legend
    legend_str = sprintf('T_%d, H_%d', i, i);
    
    % 1. Plot the actual point (This gets added to the legend)
    plot(x_dim(i), y_dim(i), 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors{i}, ...
         'MarkerEdgeColor', 'k', 'DisplayName', legend_str);
    
    % 2. Draw vertical dashed line down to the X-axis (Hidden from legend)
    plot([x_dim(i), x_dim(i)], [1e-5, y_dim(i)], '--', 'Color', colors{i}, ...
         'LineWidth', 1.2, 'HandleVisibility', 'off');
    
    % 3. Draw horizontal dashed line over to the Y-axis (Hidden from legend)
    plot([1e-3, x_dim(i)], [y_dim(i), y_dim(i)], '--', 'Color', colors{i}, ...
         'LineWidth', 1.2, 'HandleVisibility', 'off');
end

% =========================================================
% 3. Graph Formatting (Log-Log Scale)
% =========================================================
set(gca, 'XScale', 'log', 'YScale', 'log');

% Lock limits to standard Le Méhauté chart dimensions
xlim([1e-3, 1]);
ylim([1e-5, 1e-1]);
grid on;

% Add the Legend in the top-left corner (NorthWest usually has empty space here)
legend('Location', 'northwest', 'FontSize', 11);

% Professional aesthetic touches
ax = gca;
ax.LineWidth = 1.2;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.6;
ax.MinorGridLineStyle = 'none';

xlabel('Dimensionless Depth: d / (g T^2)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Dimensionless Height: H / (g T^2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Wave Theory Validity Diagram (Le Méhauté, 1976) with Evaluated Waves', 'FontSize', 14);

hold off;