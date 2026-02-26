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
T1 = t_t_next - t_t_prev;                 
T2 = t_z_down_end - t_z_down_start;       
T3 = t_z_up_end - t_z_up_start;           
T4 = t_c - t_c_prev;                      
T5 = t_c_next - t_c;                      

H1 = max_crest - min(t_prev, t_next);     
H2 = max_crest - t_prev;                  
H3 = max_crest - t_next;                  
H4 = max_crest - t_prev;                  
H5 = max_crest - t_next;                  

T = [T1, T2, T3, T4, T5];
H = [H1, H2, H3, H4, H5];

%% 5. REGULAR WAVE ANALYSIS (Wave 1 Only: T1 & H1)
d = 30;           % Water depth [m]
g = 9.81;         % Acceleration due to gravity [m/s^2]

T_eval = T(1);
H_eval = H(1);

omega = 2 * pi / T_eval;
k = fDispersion(d, T_eval); 
t_shifted = t_max_wave - t_c;

% =========================================================
% SURFACE ELEVATION CALCULATIONS (eta)
% =========================================================

% --- Linear Theory ---
eta_lin = (H_eval / 2) * cos(omega * t_shifted);
max_eta_lin = max(eta_lin);

% --- Stokes 2nd Order ---
term1 = (H_eval / 2) * cos(omega * t_shifted);
f2 = (cosh(k*d) * (2 + cosh(2*k*d))) / (16 * sinh(k*d)^3);
term2 = k * H_eval^2 * f2 * cos(2 * omega * t_shifted);
eta_s2 = term1 + term2;
max_eta_s2 = max(eta_s2);

% --- Stokes 5th Order ---
[A, B, C] = stokes_fifth_coefficients(k, d);
lambda = (H_eval * k) / 2; 
B(1,1) = 1; % Restore base linear amplitude coefficient

eta_s5 = zeros(size(t_shifted)); 
for i = 1:5      % m = power of lambda (columns)
    for j = 1:5  % n = harmonic order (rows)
        if ~isnan(B(i,j))
            eta_term = (1/k) * B(i,j) * (lambda^i) * cos(j * omega * t_shifted);
            eta_s5 = eta_s5 + eta_term;
        end
    end
end
max_eta_s5 = max(eta_s5);

% Print Errors for T1/H1 specifically
fprintf('--- CREST MATCHING ERROR (Wave 1: T1, H1) ---\n');
fprintf('Linear Error:     %.2f%%\n', 100 * abs(max_crest - max_eta_lin) / max_crest);
fprintf('Stokes 2nd Error: %.2f%%\n', 100 * abs(max_crest - max_eta_s2) / max_crest);
fprintf('Stokes 5th Error: %.2f%%\n\n', 100 * abs(max_crest - max_eta_s5) / max_crest);
crest=[max_eta_lin,max_eta_s2,max_eta_s5]
Error=[abs(max_crest - max_eta_lin),abs(max_crest - max_eta_s2),abs(max_crest - max_eta_s5)]
% =========================================================
% HORIZONTAL VELOCITY CALCULATIONS (u) AT THE CREST
% =========================================================
% Define z-vectors from the seabed (-d) up to the respective calculated crests
z_lin = linspace(-d, max_eta_lin, 50);
z_s2  = linspace(-d, max_eta_s2, 50);
z_s5  = linspace(-d, max_eta_s5, 50);

% --- Linear Theory Velocity ---
u_lin = (H_eval / 2) * omega * (cosh(k .* (z_lin + d)) ./ sinh(k * d));

% --- Stokes 2nd Order Velocity ---
u_s2_term1 = (H_eval / 2) * omega * (cosh(k .* (z_s2 + d)) ./ sinh(k * d));
u_s2_term2 = (3/4) * omega * k * (H_eval / 2)^2 * (cosh(2 * k .* (z_s2 + d)) ./ sinh(k * d)^4);
u_s2 = u_s2_term1 + u_s2_term2;

% --- Stokes 5th Order Velocity ---
% Using Matrix A (Velocity Potential coefficients)
c = omega / k;  
u_s5 = zeros(size(z_s5));

for i = 1:5      % m = power of lambda (columns)
    for j = 1:5  % n = harmonic order (rows)
        if ~isnan(A(j,i))
            % Apply derivative d/dx (pulls out a factor of n*k)
            u_term = j * k * A(j,i) * (lambda^i) * cosh(j * k .* (z_s5 + d));
            u_s5 = u_s5 + C(1) * sqrt(g/(k^3)) * u_term;
        end
    end
end

% =========================================================
% PLOTTING (Side-by-Side Layout)
% =========================================================
figure('Position', [100, 100, 1000, 500]);

% --- Subplot 1: Surface Elevation ---
subplot(1, 3, 1:2);
plot(t_shifted, eta_max_wave, 'k-', 'LineWidth', 2, 'DisplayName', 'Recorded Data');
hold on; grid on;
fprintf('eta_max=%.2fm',max_crest)
plot(t_shifted, eta_lin, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Linear Theory');
plot(t_shifted, eta_s2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Stokes 2nd');
plot(t_shifted, eta_s5, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Stokes 5th');

xlim([-5 5]); 
ylim([-5 max(max_crest)*1.2]); 
xlabel('Time [s]', 'FontWeight', 'bold'); 
ylabel('Elevation [m]', 'FontWeight', 'bold');
title(sprintf('Surface Elevation (Wave 1: T = %.2fs, H = %.2fm)', T_eval, H_eval));
legend('Location', 'northeast');
ax1 = gca; ax1.GridLineStyle = ':'; ax1.GridAlpha = 0.5;
hold off;

% --- Subplot 2: Horizontal Velocity Profile ---
subplot(1, 3, 3);
hold on; grid on;

plot(u_lin, z_lin, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Linear Theory');
plot(u_s2, z_s2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Stokes 2nd');
plot(u_s5, z_s5, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Stokes 5th');

% Dashed line for Still Water Level (SWL)
yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Velocity u [m/s]', 'FontWeight', 'bold'); 
ylabel('Elevation z [m]', 'FontWeight', 'bold');
title('Velocity Profile at Crest');
legend('Location', 'northwest');
ax2 = gca; ax2.GridLineStyle = ':'; ax2.GridAlpha = 0.5;
hold off;

[max(u_lin),max(u_s2),max(u_s5)]
