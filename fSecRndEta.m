function [eta] = fSecRndEta(x,y,t,spec)

% Second-order random wave theory, surface elevations

% eta      - [m] free surface
% x, y     - [m] spatial coordinates
% t        - [s] time instances
% spec     - [] wave spectrum, after being processed by fSpecImport

% x, y, t should be one of the following:
% 1) All scalars.
% 2) One of them is an array, others are scalars
% 3) All of them are arrays, with identical dimensions.

% Li Ma (2013), lm808@ic.ac.uk, Imperial College London


Nc = length(spec.f);

eta1 = zeros(size(x+y+t));
eta2 = zeros(size(eta1));

%% SpecPrep
R = spec.k .* tanh(spec.k*spec.d);
B1 = 0.25 * sum(spec.a.^2 .* (2*R - (spec.k.^2+R.^2)./R));

%% Double summation
for i = 1:Nc-1
    dgt = length([num2str(i-1),' /',num2str(Nc)]);
    if i>1
        fprintf(1,repmat('\b',1,dgt));
    end
    fprintf(1,' %u/%u',i,Nc)

    Diip = (12*R(i)*(spec.k(i)^2-R(i)^2))./(4*R(i)-2*spec.k(i)*tanh(2*spec.k(i)*spec.d));
    B2 = 0.25 * spec.a(i)^2 * (Diip/R(i) + 2*R(i));
    
    THi = (spec.k(i) * (x*cos(spec.th(i)) + y*sin(spec.th(i))) - spec.omega(i)*t) + spec.ph(i);
    eta1 = eta1 + spec.a(i) * cos(THi);
    eta2 = eta2 + B2 * cos(2*THi);
    
    for j = (i+1):Nc
        gamma = cos(spec.th(i)-spec.th(j));
        kp = sqrt(spec.k(i)^2 + spec.k(j)^2 + 2*spec.k(i)*spec.k(j)*gamma);
        km = sqrt(spec.k(i)^2 + spec.k(j)^2 - 2*spec.k(i)*spec.k(j)*gamma);
        Dijm = (((sqrt(R(i))-sqrt(R(j)))*(sqrt(R(j))*(spec.k(i)^2-R(i)^2)-sqrt(R(i))*(spec.k(j)^2-R(j)^2)))...
               +(2*((sqrt(R(i))-sqrt(R(j)))^2)*(spec.k(i)*spec.k(j)*gamma+R(i)*R(j))))...
               /(((sqrt(R(i))-sqrt(R(j)))^2)-km*tanh(km*spec.d));
        Dijp = (((sqrt(R(i))+sqrt(R(j)))*(sqrt(R(j))*(spec.k(i)^2-R(i)^2)+sqrt(R(i))*(spec.k(j)^2-R(j)^2)))...
               +(2*((sqrt(R(i))+sqrt(R(j)))^2)*(spec.k(i)*spec.k(j)*gamma-R(i)*R(j))))...
               /(((sqrt(R(i))+sqrt(R(j)))^2)-kp*tanh(kp*spec.d));
        
        B3 = 0.5*spec.a(i)*spec.a(j)* ((Dijm-(spec.k(i)*spec.k(j)*gamma + R(i)*R(j)))/sqrt(R(i)*R(j))+R(i)+R(j));
        B4 = 0.5*spec.a(i)*spec.a(j)* ((Dijp-(spec.k(i)*spec.k(j)*gamma - R(i)*R(j)))/sqrt(R(i)*R(j))+R(i)+R(j));
        
        THj = (spec.k(j) * (x*cos(spec.th(j)) + y*sin(spec.th(j))) - spec.omega(j)*t) + spec.ph(j);
        eta2 = eta2 + B3 * cos(THi-THj) + B4 * cos(THi+THj);
        
        if isnan(gamma)
            error('')
        end
    end
end
i = Nc;
Diip = (12*R(i)*(spec.k(i)^2-R(i)^2))./(4*R(i)-2*spec.k(i)*tanh(2*spec.k(i)*spec.d));
B2 = 0.25 * spec.a(i)^2 * (Diip/R(i) + 2*R(i));
THi = (spec.k(i) * (x*cos(spec.th(i)) + y*sin(spec.th(i))) - spec.omega(i)*t) + spec.ph(i);

eta1 = eta1 + spec.a(i) * cos(THi);
eta2 = eta2 + B2 * cos(2*THi);

dgt = length([num2str(i-1),' /',num2str(Nc)]);
fprintf(1,repmat('\b',1,dgt));

%%
eta2 = eta2 + B1;
eta = eta1 +  eta2;
