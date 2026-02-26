function [u,v,w] = fSecRndVel(x,y,z,t,spec)

% Second-order random wave theory, velocities

% u, v, w  - [m/s] velocity in x, y, z directions
% x, y, z  - [m]   spatial coordinates
% t        - [s]   time instances
% spec     - []    wave spectrum, after being processed by fSpecImport

% x, y, z, t should be one of the following:
% 1) All scalars.
% 2) One of them is an array, others are scalars
% 3) All of them are arrays, with identical dimensions.

% Li Ma (2013), lm808@ic.ac.uk, Imperial College London


Nc = length(spec.f);
d = spec.d;

u1 = zeros(size(x+y+z+t));
u2 = u1;

v1 = u1;
v2 = u1;
 
w1 = u1;
w2 = u1;

g = 9.81;
R = spec.k .* tanh(spec.k*d);

%% Double summation
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
    
    u1 = u1 + A0 * spec.k(i)*cos(spec.th(i)) * cosh(spec.k(i)*(z+d)) .* cos(THi);
    u2 = u2 + A1 * 2*spec.k(i)*cos(spec.th(i)) * cosh(2*spec.k(i)*(z+d)) .* cos(2*THi);

    v1 = v1 + A0 * spec.k(i)*sin(spec.th(i)) * cosh(spec.k(i)*(z+d)) .* cos(THi);
    v2 = v2 + A1 * 2*spec.k(i)*sin(spec.th(i)) * cosh(2*spec.k(i)*(z+d)) .* cos(2*THi);

    w1 = w1 + A0 * spec.k(i) *sinh(spec.k(i)*(z+d)) .* sin(THi);
    w2 = w2 + A1 * 2*spec.k(i) * sinh(2*spec.k(i)*(z+d)) .* sin(2*THi);
    
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
        
        u2 = u2...
             + A2 * (spec.k(i)*cos(spec.th(i))-spec.k(j)*cos(spec.th(j))) * cosh(km*(z+d)) .* cos(THi-THj)...
             + A3 * (spec.k(i)*cos(spec.th(i))+spec.k(j)*cos(spec.th(j))) * cosh(kp*(z+d)) .* cos(THi+THj);
         
        v2 = v2...
             + A2 * (spec.k(i)*sin(spec.th(i))-spec.k(j)*sin(spec.th(j))) * cosh(km*(z+d)) .* cos(THi-THj)...
             + A3 * (spec.k(i)*sin(spec.th(i))+spec.k(j)*sin(spec.th(j))) * cosh(kp*(z+d)) .* cos(THi+THj);
         
        w2 = w2...
             + A2 * km * sinh(km*(z+d)) .* sin(THi-THj)...
             + A3 * kp * sinh(kp*(z+d)) .* sin(THi+THj);
    end
end

i = Nc;

Diip = (12*R(i)*(spec.k(i)^2-R(i)^2))./(4*R(i)-2*spec.k(i)*tanh(2*spec.k(i)*d));
A0 = spec.a(i) * g / spec.omega(i) / cosh(spec.k(i)*d);
A1 = 0.125 * (spec.a(i)^2) *g*g * Diip / (spec.omega(i)^3) / cosh(2*spec.k(i)*d);
    
THi = (spec.k(i) * (x*cos(spec.th(i)) + y*sin(spec.th(i))) - spec.omega(i)*t) + spec.ph(i);

u1 = u1 + A0 * spec.k(i)*cos(spec.th(i)) * cosh(spec.k(i)*(z+d)) .* cos(THi);
u2 = u2 + A1 * 2*spec.k(i)*cos(spec.th(i)) * cosh(2*spec.k(i)*(z+d)) .* cos(2*THi);

v1 = v1 + A0 * spec.k(i)*sin(spec.th(i)) * cosh(spec.k(i)*(z+d)) .* cos(THi);
v2 = v2 + A1 * 2*spec.k(i)*sin(spec.th(i)) * cosh(2*spec.k(i)*(z+d)) .* cos(2*THi);

w1 = w1 + A0 * spec.k(i) *sinh(spec.k(i)*(z+d)) .* sin(THi);
w2 = w2 + A1 * 2*spec.k(i) * sinh(2*spec.k(i)*(z+d)) .* sin(2*THi);
    
dgt = length([num2str(i-1),' /',num2str(Nc)]);
fprintf(1,repmat('\b',1,dgt));

%%
u = u1+u2;
v = v1+v2;
w = w1+w2;