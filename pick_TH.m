clear; clc; close all;

%% 1. Load the data
load('Data_Fig1.mat'); 

%% 2. Identify Crests, Troughs, and Zero-Crossings
[crests, t_crests] = findpeaks(eta_max_wave, t_max_wave);
[troughs_inv, t_troughs] = findpeaks(-eta_max_wave, t_max_wave);
troughs = -troughs_inv;

idx_zero = find(eta_max_wave(1:end-1) .* eta_max_wave(2:end) < 0);
t_zero = t_max_wave(idx_zero) - eta_max_wave(idx_zero) .* (t_max_wave(idx_zero+1) - t_max_wave(idx_zero)) ./ (eta_max_wave(idx_zero+1) - eta_max_wave(idx_zero));

is_up_crossing = eta_max_wave(idx_zero) < 0; 
t_z_up = t_zero(is_up_crossing);
t_z_down = t_zero(~is_up_crossing);

%% 3. Extract Specific Points Around the Main Wave
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

%% 4. Evaluate T and H
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

%% 5. Linear Regular Theory
t_shifted = t_max_wave - t_c;

figure('Position', [150, 150, 800, 600]);
subplot(1,3,1)
plot(t_shifted, eta_max_wave, 'k-', 'LineWidth', 2, 'DisplayName', 'Recorded Data');
hold on; grid on;

colors = {'b', 'r', [0 0.5 0], 'm', [0.85 0.33 0.1]};
err_pct_linear = zeros(1,5);

for i = 1:5
    omega = 2 * pi / T(i);
    eta_linear = (H(i) / 2) * cos(omega * t_shifted);
    
    % Calculate Crest Error
    eta_max_theory = max(eta_linear); 
    crest_err = abs(max_crest - eta_max_theory);
    err_pct_linear(i) = 100 * crest_err / max_crest;
    
    plot_label = sprintf('Linear Theory (T_%d, H_%d)', i, i);
    plot(t_shifted, eta_linear, '-', 'Color', colors{i}, 'LineWidth', 1.5, 'DisplayName', plot_label);
end

% Print the best match
[best_err_lin, best_idx_lin] = min(err_pct_linear);
fprintf('--- LINEAR THEORY ---\n');
fprintf('Best match: Wave %d (Error: %.2f%%)\n\n', best_idx_lin, best_err_lin);

xlim([-5 5]); ylim([-5 8]); 
xlabel('Time [s]'); ylabel('Elevation [m]');
title('Surface Elevation (Linear Theory)');
legend('Location', 'northeast');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5; hold off;

%% 6. Second-Order Regular Theory
d = 30;           % Water depth [m]
g = 9.81;         % Acceleration due to gravity [m/s^2]


subplot(1,3,2)
plot(t_shifted, eta_max_wave, 'k-', 'LineWidth', 2, 'DisplayName', 'Recorded Data');
hold on; grid on;

err_pct_stokes2 = zeros(1,5);

for i = 1:5
    T_current = T(i);
    H_cur = H(i);
    omega = 2 * pi / T_current;
    
    % Fixed inputs: d first, T second
    k = fDispersion(d, T_current); 
    
    term1 = (H_cur / 2) * cos(omega * t_shifted);
    f2 = (cosh(k*d) * (2 + cosh(2*k*d))) / (16 * sinh(k*d)^3);
    term2 = k * H_cur^2 * f2 * cos(2 * omega * t_shifted);
    eta_stokes2 = term1 + term2;
    
    % Calculate Crest Error
    eta_max_theory = max(eta_stokes2); 
    crest_err = abs(max_crest - eta_max_theory);
    err_pct_stokes2(i) = 100 * crest_err / max_crest;
    
    plot_label = sprintf('Stokes 2nd (T_%d, H_%d)', i, i);
    plot(t_shifted, eta_stokes2, '-', 'Color', colors{i}, 'LineWidth', 1.5, 'DisplayName', plot_label);
end

% Print the best match
[best_err_s2, best_idx_s2] = min(err_pct_stokes2);
fprintf('--- STOKES 2ND ORDER ---\n');
fprintf('Best match: Wave %d (Error: %.2f%%)\n\n', best_idx_s2, best_err_s2);

xlim([-5 5]); ylim([-5 8]); 
xlabel('Time [s]'); ylabel('Elevation [m]');
title('Surface Elevation (Stokes 2nd Order)');
legend('Location', 'northeast');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5; hold off;

%% 7. Fifth Order Regular Theory\
subplot(1,3,3)

plot(t_shifted, eta_max_wave, 'k-', 'LineWidth', 2, 'DisplayName', 'Recorded Data');
hold on; grid on;

err_pct_stokes5 = zeros(1,5);

for w = 1:5  
    T_current = T(w);
    H_cur = H(w);
    
    k = fDispersion(d, T_current); 
    omega = 2 * pi / T_current;
    
    [A, B, C] = stokes_fifth_coefficients(k, d);
    
    lambda = (H_cur * k) / 2; 
    
    % Restore the missing base linear term
    B(1,1) = 1; 
    
    eta_stokes5 = zeros(size(t_shifted)); 
    
    for i = 1:5
        for j = 1:5
            if ~isnan(B(i,j))
                eta_term = (1/k) * B(i,j) * (lambda^i) * cos(j * omega * t_shifted);
                eta_stokes5 = eta_stokes5 + eta_term;
            end
        end
    end
    
    % Calculate Crest Error
    eta_max_theory = max(eta_stokes5); 
    crest_err = abs(max_crest - eta_max_theory);
    err_pct_stokes5(w) = 100 * crest_err / max_crest;
    
    plot_label = sprintf('Stokes 5th (T_%d, H_%d)', w, w);
    plot(t_shifted, eta_stokes5, '-', 'Color', colors{w}, 'LineWidth', 1.5, 'DisplayName', plot_label);
end

% Print the best match
[best_err_s5, best_idx_s5] = min(err_pct_stokes5);
fprintf('--- STOKES 5TH ORDER ---\n');
fprintf('Best match: Wave %d (Error: %.2f%%)\n\n', best_idx_s5, best_err_s5);

xlim([-5 5]); ylim([-5 8]); 
xlabel('Time [s]'); ylabel('Elevation [m]');
title('Surface Elevation (Stokes 5th Order)');
legend('Location', 'northeast');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5; hold off;