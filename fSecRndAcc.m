function [dudt,dvdt,dwdt] = fSecRndAcc(x,y,z,t,spec)

% Second-order random wave theory, accelerations

% dudt, dvdt, dwdt  - [m/s^2] acceleration in x, y, z directions
% x, y, z           - [m]     spatial coordinates
% t                 - [s]     time instances
% spec              - []      wave spectrum, after being processed by fSpecImport

% x, y, z, t should be one of the following:
% 1) All scalars.
% 2) One of them is an array, others are scalars
% 3) All of them are arrays, with identical dimensions.

% Li Ma (2013), lm808@ic.ac.uk, Imperial College London


Nc = length(spec.f);
d = spec.d;

dudt1 = zeros(size(x+y+z+t));
dudt2 = dudt1;

dvdt1 = dudt1;
dvdt2 = dudt1;
 
dwdt1 = dudt1;
dwdt2 = dudt1;

%%
g = 9.81;
R = spec.k .* tanh(spec.k*d);

for i = 1:(Nc-1)

    dgt = length([num2str(i-1),' /',num2str(Nc)]);
    if i>1
        fprintf(1,repmat('\b',1,dgt));
    end
    fprintf(1,' %u/%u',i,Nc)
    
    Diip = (12*R(i)*(spec.k(i)^2-R(i)^2))./(4*R(i)-2*spec.k(i)*tanh(2*spec.k(i)*d));
    A0 = spec.a(i) * g / spec.omega(i) / cosh(spec.k(i)*d);
    A1 = 0.125 * (spec.a(i)^2) *g*g * Diip / (spec.omega(i)^3) / cosh(2*spec.k(i)*d);

    THi = (spec.k(i) * (x*cos(spec.th(i)) + y*sin(spec.th(i))) - spec.omega(i)*t) + spec.ph(i);
    
    dudt1 = dudt1 + A0 * spec.k(i)*cos(spec.th(i)) * spec.omega(i) * cosh(spec.k(i)*(z+d)) .* sin(THi);
    dudt2 = dudt2 + A1 * 2*spec.k(i)*cos(spec.th(i)) * 2*spec.omega(i) * cosh(2*spec.k(i)*(z+d)) .* sin(2*THi);
    
    dvdt1 = dvdt1 + A0 * spec.k(i)*sin(spec.th(i)) * spec.omega(i) * cosh(spec.k(i)*(z+d)) .* sin(THi);
    dvdt2 = dvdt2 + A1 * 2*spec.k(i)*sin(spec.th(i)) * 2*spec.omega(i) * cosh(2*spec.k(i)*(z+d)) .* sin(2*THi);
    
    dwdt1 = dwdt1 + A0 * spec.k(i) * (-spec.omega(i)) * sinh(spec.k(i)*(z+d)) .* cos(THi);
    dwdt2 = dwdt2 + A1 * 2*spec.k(i) * (-2*spec.omega(i)) * sinh(2*spec.k(i)*(z+d)) .* cos(2*THi);
    
    for j = (i+1):Nc
        
        gamma = cos(spec.th(i)-spec.th(j));
        kp = sqrt(spec.k(i)^2 + spec.k(j)^2 + 2*spec.k(i)*spec.k(j)*gamma);
        km = sqrt(spec.k(i)^2 + spec.k(j)^2 - 2*spec.k(i)*spec.k(j)*gamma);
        Dijm = (((sqrt(R(i))-sqrt(R(j)))*(sqrt(R(j))*(spec.k(i)^2-R(i)^2)-sqrt(R(i))*(spec.k(j)^2-R(j)^2)))...
               +(2*((sqrt(R(i))-sqrt(R(j)))^2)*(spec.k(i)*spec.k(j)*gamma+R(i)*R(j))))...
               /(((sqrt(R(i))-sqrt(R(j)))^2)-km*tanh(km*d));
        Dijp = (((sqrt(R(i))+sqrt(R(j)))*(sqrt(R(j))*(spec.k(i)^2-R(i)^2)+sqrt(R(i))*(spec.k(j)^2-R(j)^2)))...
               +(2*((sqrt(R(i))+sqrt(R(j)))^2)*(spec.k(i)*spec.k(j)*gamma-R(i)*R(j))))...
               /(((sqrt(R(i))+sqrt(R(j)))^2)-kp*tanh(kp*d));

        A2 = 0.5 * spec.a(i)*spec.a(j)*g*g*Dijm / (spec.omega(i)*spec.omega(j)*(spec.omega(i)-spec.omega(j))*cosh(km*d));
        A3 = 0.5 * spec.a(i)*spec.a(j)*g*g*Dijp / (spec.omega(i)*spec.omega(j)*(spec.omega(i)+spec.omega(j))*cosh(kp*d));
        
        if spec.omega(i)==spec.omega(j)
            A2 = 0;
        end
        
        THj = (spec.k(j) * (x*cos(spec.th(j)) + y*sin(spec.th(j))) - spec.omega(j)*t) + spec.ph(j);
        
        dudt2 = dudt2...
             + A2 * (spec.k(i)*cos(spec.th(i))-spec.k(j)*cos(spec.th(j))) * (spec.omega(i)-spec.omega(j)) * cosh(km*(z+d)) .* sin(THi-THj)...
             + A3 * (spec.k(i)*cos(spec.th(i))+spec.k(j)*cos(spec.th(j))) * (spec.omega(i)+spec.omega(j)) * cosh(kp*(z+d)) .* sin(THi+THj);
         
        dvdt2 = dvdt2...
             + A2 * (spec.k(i)*sin(spec.th(i))-spec.k(j)*sin(spec.th(j))) * (spec.omega(i)-spec.omega(j)) * cosh(km*(z+d)) .* sin(THi-THj)...
             + A3 * (spec.k(i)*sin(spec.th(i))+spec.k(j)*sin(spec.th(j))) * (spec.omega(i)+spec.omega(j)) * cosh(kp*(z+d)) .* sin(THi+THj);
         
        dwdt2 = dwdt2...
             + A2 * km * (-spec.omega(i)+spec.omega(j)) * sinh(km*(z+d)) .* cos(THi-THj)...
             + A3 * kp * (-spec.omega(i)-spec.omega(j)) * sinh(kp*(z+d)) .* cos(THi+THj);
    end
end

i = Nc;

Diip = (12*R(i)*(spec.k(i)^2-R(i)^2))./(4*R(i)-2*spec.k(i)*tanh(2*spec.k(i)*d));
A0 = spec.a(i) * g / spec.omega(i) / cosh(spec.k(i)*d);
A1 = 0.125 * (spec.a(i)^2) *g*g * Diip / (spec.omega(i)^3) / cosh(2*spec.k(i)*d);

THi = (spec.k(i) * (x*cos(spec.th(i)) + y*sin(spec.th(i))) - spec.omega(i)*t) + spec.ph(i);

dudt1 = dudt1 + A0 * spec.k(i)*cos(spec.th(i)) * spec.omega(i) * cosh(spec.k(i)*(z+d)) .* sin(THi);
dudt2 = dudt2 + A1 * 2*spec.k(i)*cos(spec.th(i)) * 2*spec.omega(i) * cosh(2*spec.k(i)*(z+d)) .* sin(2*THi);

dvdt1 = dvdt1 + A0 * spec.k(i)*sin(spec.th(i)) * spec.omega(i) * cosh(spec.k(i)*(z+d)) .* sin(THi);
dvdt2 = dvdt2 + A1 * 2*spec.k(i)*sin(spec.th(i)) * 2*spec.omega(i) * cosh(2*spec.k(i)*(z+d)) .* sin(2*THi);

dwdt1 = dwdt1 + A0 * spec.k(i) * (-spec.omega(i)) * sinh(spec.k(i)*(z+d)) .* cos(THi);
dwdt2 = dwdt2 + A1 * 2*spec.k(i) * (-2*spec.omega(i)) * sinh(2*spec.k(i)*(z+d)) .* cos(2*THi);
    
dgt = length([num2str(i-1),' /',num2str(Nc)]);
fprintf(1,repmat('\b',1,dgt));

%%
dudt = dudt1+dudt2;
dvdt = dvdt1+dvdt2;
dwdt = dwdt1+dwdt2;