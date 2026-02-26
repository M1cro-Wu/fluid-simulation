% Solves the wave dispersion equation
% h     - depth [m]
% T     - wave period [s] 
% k     - wavenumber [1/m]
% JS, March 2010

function k = fDispersion(h,T)

k_all = zeros(size(T));

for i = 1:length(T)
    g = 9.81;
    omega = 2*pi/T(i);

    p = omega^2*h/g;
    q = (tanh(p^0.75))^(-2/3);
    k0 = omega^2*q/g;

    [k,F] = fzero(@(k) omega^2 -g*k*tanh(k*h),k0);
    
    k_all(i) = k;
end

k = k_all;