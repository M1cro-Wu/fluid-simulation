clear; clc; close all;

%% 1. Load Data
load('Data_Fig2.mat'); % Random Sea Spectrum
load('Data_Fig1.mat'); % Recorded Freak Wave
d = 30; 
g = 9.81; 

% --- FULL SPECTRUM VARIABLES ---
f_full = f(2:end);
phi_full = phase(2:end); 
a_full = amp(2:end);     

N_full = length(f_full);
omega_full = 2 * pi * f_full;
T_full = 1./f_full;

% Calculate wavenumbers for FULL spectrum
k_full = zeros(N_full, 1);
for n = 1:N_full
    k_full(n) = fDispersion(d, T_full(n)); 
end

%% 2. Generate 1st Order Elevation (Full Spectrum)
t = linspace(-60, 60, 2400); % Fixed time vector
eta1 = zeros(size(t));

for n = 1:N_full
    eta1 = eta1 + a_full(n) * cos(omega_full(n) * t + phi_full(n));
end

%% 3. Velocity Profiles (Full Spectrum)
[max_eta1, ~] = max(eta1);
z_vec = linspace(-d, max_eta1, 100);
u_wheeler = zeros(size(z_vec));
u_linear = zeros(size(z_vec)); 

% Wheeler mapped coordinate
z_prime = (d * (z_vec + d) / (max_eta1 + d)) - d;

for n = 1:N_full
    % Wheeler Stretching
    ratio_wheeler = cosh(k_full(n) * (z_prime + d)) / sinh(k_full(n) * d);
    if k_full(n)*d > 50; ratio_wheeler = exp(k_full(n)*z_prime); end
    u_wheeler = u_wheeler + a_full(n) * omega_full(n) * ratio_wheeler * cos(phi_full(n));
    
    % Pure Linear (No Stretch)
    ratio_linear = cosh(k_full(n) * (z_vec + d)) / sinh(k_full(n) * d);
    if k_full(n)*d > 50; ratio_linear = exp(k_full(n)*z_vec); end
    u_linear = u_linear + a_full(n) * omega_full(n) * ratio_linear * cos(phi_full(n));
end

%% 4. 2nd Order Calculation (Filtered)

% amp_threshold = 0.01 * max(a_full); 
keep_idx = f >= 0.05 & f <= 0.3;

% Create FILTERED subsets
f_filt = f_full(keep_idx);
a_filt = a_full(keep_idx);
phi_filt = phi_full(keep_idx);
k_filt = k_full(keep_idx);         
omega_filt = omega_full(keep_idx); 

% Build the 'spec' struct for the function
spec.f = f_filt;
spec.a = a_filt;
spec.ph = phi_filt;
spec.omega = omega_filt;
spec.k = k_filt;
spec.d = d;
spec.th = zeros(size(f_filt));

% Calculate 2nd Order Total (Filtered)
x = 0; y = 0;
eta_2nd_filtered_total = fSecRndEta(x, y, t, spec);



% C. Calculate 2nd Order Velocity Profile (Using fSecRndVel)
% Note: fSecRndVel returns (u, v, w). We only need u (horizontal).
[u_2nd, ~, ~] = fSecRndVel(0, 0, z_vec, 0, spec);

% Isolate correction term: (Total Filtered) - (Linear Filtered)
eta_1st_filtered_part = zeros(size(t));
for n = 1:length(f_filt)
    eta_1st_filtered_part = eta_1st_filtered_part + a_filt(n) * cos(omega_filt(n) * t + phi_filt(n));
end
eta_2nd_correction = eta_2nd_filtered_total - eta_1st_filtered_part;

% Add correction to FULL linear signal
eta_final_hybrid = eta1 + eta_2nd_correction;

%% 5. Peak Alignment
% Align Recorded Data
[max_rec, idx_rec] = max(eta_max_wave);
t_rec_shifted = t_max_wave - t_max_wave(idx_rec);

% Align 1st Order
[~, idx_1st] = max(eta1);
t_sim_shifted = t - t(idx_1st);

% Align 2nd Order Hybrid
[max_2nd, idx_2nd] = max(eta_final_hybrid);
t_2nd_shifted = t - t(idx_2nd);

%% 6. PLOTTING (Sample Code Layout)
figure('Position', [100, 100, 1200, 500]);

% --- Subplot 1: Surface Elevation (Spans 2 columns) ---
subplot(1, 3, 1:2); 
plot(t_rec_shifted, eta_max_wave, 'k-', 'LineWidth', 2, 'DisplayName', 'Recorded Data');
hold on; grid on;
plot(t_sim_shifted, eta1, 'b--', 'LineWidth', 1.5, 'DisplayName', '1st Order (Linear)');
plot(t_2nd_shifted, eta_final_hybrid, 'r-', 'LineWidth', 1.5, 'DisplayName', '2nd Order (Non-linear)');

xlim([-10 10]);
ylim([min(eta_max_wave)*1.2, max(eta_max_wave)*1.2]);
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Elevation \eta [m]', 'FontWeight', 'bold');
title('Random Wave Elevation Comparison');
legend('Location', 'northeast');

% --- Subplot 2: Velocity Profile (Spans 1 column) ---
subplot(1, 3, 3);
plot(u_wheeler, z_vec, 'b-', 'LineWidth', 2, 'DisplayName', 'Wheeler Stretch');
hold on; grid on;
plot(u_linear, z_vec, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Linear (No Stretch)');
plot(u_2nd, z_vec, 'r-', 'LineWidth', 2, 'DisplayName', '2nd Order (Total)');


yline(0, 'k-', 'HandleVisibility', 'off'); % MWL
yline(max_eta1, 'k:', 'HandleVisibility', 'off'); % Crest Height
xlim([0 10]);

xlabel('Velocity u [m/s]', 'FontWeight', 'bold');
ylabel('Elevation z [m]', 'FontWeight', 'bold');
title('Horizontal Velocity Profile');
legend('Location', 'southeast');

%% 7. Error Output
fprintf('\n--- CREST ERROR ---\n');
fprintf('1st Order Error: %.2f%%\n', 100 * abs(max_rec - max(eta1)) / max_rec);
fprintf('2nd Order Error: %.2f%%\n', 100 * abs(max_rec - max(eta_final_hybrid)) / max_rec);
crest=[max(eta1),max(eta_final_hybrid)]
Error=abs(max_rec-crest)
max(u_linear),max(u_wheeler),max(u_2nd)