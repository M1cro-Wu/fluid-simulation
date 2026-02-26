function [k, c, L] = stokes5th_dispersion(T, H, d, g)
% STOKES5TH_DISPERSION  Solve for wave number k using Fenton's (1985)
% fifth-order Stokes wave theory dispersion relation.
%
% Reference: Fenton, J.D. (1985). "A Fifth-Order Stokes Theory for 
%            Steady Waves." Journal of Waterway, Port, Coastal and 
%            Ocean Engineering, 111(2), 216-234.
%
% INPUTS:
%   T  - Wave period [s]
%   H  - Wave height [m]
%   d  - Water depth [m]
%   g  - Gravitational acceleration [m/s^2] (default: 9.81)
%
% OUTPUTS:
%   k  - Wave number [rad/m]
%   c  - Phase speed [m/s]
%   L  - Wavelength [m]
%
% USAGE:
%   [k, c, L] = stokes5th_dispersion(10, 2, 20)
%   [k, c, L] = stokes5th_dispersion(10, 2, 20, 9.81)

    if nargin < 4
        g = 9.81;
    end

    omega = 2 * pi / T;  % Angular frequency

    %% ===================================================================
    %  STEP 1: Initial estimate from linear dispersion relation
    %          omega^2 = g*k*tanh(k*d)
    %  Solved iteratively with Newton-Raphson
    % ====================================================================
    k0 = omega^2 / g;  % Deep-water initial guess
    for i = 1:100
        f  = omega^2 - g * k0 * tanh(k0 * d);
        df = -g * (tanh(k0 * d) + k0 * d * sech(k0 * d)^2);
        dk = -f / df;
        k0 = k0 + dk;
        if abs(dk / k0) < 1e-12
            break
        end
    end
    fprintf('Linear dispersion k0 = %.6f rad/m (converged in %d iterations)\n', k0, i);

    %% ===================================================================
    %  STEP 2: Compute Fenton (1985) fifth-order coefficients
    %  Using Stokes' first definition of wave celerity (zero mean current)
    % ====================================================================
    k = k0;  % Start with linear estimate

    for iter = 1:200
        S  = 1 / cosh(2 * k * d);  % S = sech(2kd)
        S2 = S^2;
        S4 = S^4;

        % Shorthand
        kd = k * d;
        C  = cosh(kd);
        Sh = sinh(kd);
        t  = tanh(kd);

        % -----------------------------------------------------------------
        %  Coefficients for the dispersion relation (Fenton 1985, Table 1)
        %  c^2 = (g/k)*tanh(kd) * [1 + eps^2*C1 + eps^4*C2]
        %  where eps = k*H/2
        % -----------------------------------------------------------------

        % C1 coefficient (second-order correction)
        C1 = (2 + 7*S2) / (4*(1 - S)^2);

        % C2 coefficient (fourth-order correction)
        C2 = (4 + 32*S - 116*S2 - 400*S^3 - 71*S4 + 146*S^5) ...
             / (32 * (1 - S)^5);

        % Wave steepness parameter
        eps = k * H / 2;

        % Fifth-order dispersion relation
        % omega^2 = g*k*tanh(kd) * [1 + eps^2*C1 + eps^4*C2]
        c2_theory = (g / k) * tanh(kd) * (1 + eps^2 * C1 + eps^4 * C2);
        c_theory  = sqrt(c2_theory);

        % The wave speed from period and wavenumber must satisfy c = omega/k
        c_target = omega / k;

        % Residual: F(k) = c_theory - c_target = 0
        F = c_theory - c_target;

        % Numerical derivative dF/dk (central difference)
        dk_num = k * 1e-7;
        k_plus  = k + dk_num;
        k_minus = k - dk_num;

        % Evaluate at k + dk
        eps_p  = k_plus * H / 2;
        S_p    = 1 / cosh(2 * k_plus * d);
        C1_p   = (2 + 7*S_p^2) / (4*(1 - S_p)^2);
        C2_p   = (4 + 32*S_p - 116*S_p^2 - 400*S_p^3 - 71*S_p^4 + 146*S_p^5) ...
                 / (32 * (1 - S_p)^5);
        c_p    = sqrt((g / k_plus) * tanh(k_plus * d) * ...
                 (1 + eps_p^2 * C1_p + eps_p^4 * C2_p));
        F_plus = c_p - omega / k_plus;

        % Evaluate at k - dk
        eps_m  = k_minus * H / 2;
        S_m    = 1 / cosh(2 * k_minus * d);
        C1_m   = (2 + 7*S_m^2) / (4*(1 - S_m)^2);
        C2_m   = (4 + 32*S_m - 116*S_m^2 - 400*S_m^3 - 71*S_m^4 + 146*S_m^5) ...
                 / (32 * (1 - S_m)^5);
        c_m    = sqrt((g / k_minus) * tanh(k_minus * d) * ...
                 (1 + eps_m^2 * C1_m + eps_m^4 * C2_m));
        F_minus = c_m - omega / k_minus;

        dFdk = (F_plus - F_minus) / (2 * dk_num);

        % Newton-Raphson update
        dk = -F / dFdk;
        k  = k + dk;

        if abs(dk / k) < 1e-12
            fprintf('Stokes 5th-order converged in %d iterations\n', iter);
            break
        end
    end

    if iter == 200
        warning('Stokes 5th-order did not converge after 200 iterations.');
    end

    %% ===================================================================
    %  STEP 3: Compute final outputs
    % ====================================================================
    c = omega / k;
    L = 2 * pi / k;

    %% ===================================================================
    %  STEP 4: Print summary
    % ====================================================================
    fprintf('\n========== Stokes 5th-Order Dispersion Results ==========\n');
    fprintf('  Input:  T = %.2f s,  H = %.2f m,  d = %.2f m\n', T, H, d);
    fprintf('  ---------------------------------------------------------\n');
    fprintf('  Wave number      k = %.6f rad/m\n', k);
    fprintf('  Wavelength       L = %.4f m\n', L);
    fprintf('  Phase speed      c = %.4f m/s\n', c);
    fprintf('  Steepness   kH/2   = %.6f  (eps)\n', k * H / 2);
    fprintf('  Relative depth kd  = %.4f\n', k * d);
    fprintf('  Ursell number U    = %.4f\n', (H * L^2) / d^3);
    fprintf('=========================================================\n');

    % Validity check
    eps_final = k * H / 2;
    kd_final  = k * d;
    if eps_final > 0.5
        warning('Wave steepness eps = %.3f is large. Stokes theory may be unreliable.', eps_final);
    end
    if kd_final < 1.0
        warning('kd = %.3f is small (shallow water). Consider cnoidal or stream function theory.', kd_final);
    end
end
