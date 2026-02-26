 % Sample codes for the application of the second-order random wave theory
% Year 2 Fluid Mechanics, Department of Civil and Environmental Engineering
% Li Ma (2017), lm808@ic.ac.uk, Imperial College London

% The following is just an example.
% You should carefully read the comments at the beginning of each function.

clear; clc; close all

%% Define constants

% water depth
d = 1.5; 

% spatial coordinates
x = 0;
y = 0;

% time vector
t = -2:0.05:2;

%% Load the wave spectrum

load sample_spectrum.mat 

% Look into the above file and figure out what kind of input format you
% need for the wave spectrum. 
% Hint: check out the MATLAB data type, struct.

spec = fSpecImport(spec, d, 1, 0); % do not worry about the last 2 inputs, keep them as 1 and 0.

%% Calculate surface
[eta] = fSecRndEta(x, y, t, spec);

%% Calculate velocity beneath the maximum surface elevation

% For a focused wave (as is the case now), it will be at x=0, y=0, t=0.
% In other cases, you may need to first find the space-time coordinates of 
% the maximum surface elevation.

% Construct a vector that extends from the bed to the surface at one 
% spatial point in plan
z = linspace(-d, max(eta), 20);

% Compute the velocities
[u, v, w] = fSecRndVel(0, 0, z, 0, spec);

% Accelerations can be calculated the same way using the function fSecRndAcc.m
% Think where the maximum acceleration may occur. Is it still where the
% maximum surface elevation is?

%% Plot the results

figure
subplot(1, 3, 1:2)
plot(t, eta)
xlabel ('Time [s]')
ylabel ('Elevation [m]')
title ('Surface elevation')
grid on; grid minor
subplot(1, 3, 3)
plot(u, z)
xlabel ('Velocity [m/s]')
ylabel ('Elevation [m]')
title ('Water particle velocity')
grid on; grid minor

