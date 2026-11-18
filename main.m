%{
===========================================================================
 File Name:    Prior_Project4.m
 Author:       Kobe Prior
 Date:         November 17, 2026
 Description:  The purpose of this program is to solve for the reflection
               and transmission of electric fields through an arbitrary 
               multi-layer medium using the 1D FDFD method.

 Notes:        assuming exp(+jwt)
===========================================================================
%}
clear; clc;

%% Input Parameters
run_layer_optimization = true;
% Domain size and discretization
domainLength  = 1.0;                  % meters
numGridPoints = 20000;
x  = linspace(0, domainLength, numGridPoints);

% Material properties for each layer
% Format: [Material 1, Material 2, ..., Material N]
eps_r_list   = [4.0, 6.0];
mu_r_list    = [4.0, 6.0];
sigma_e_list = [0.0, 0.0];
sigma_m_list = [0.0, 0.0];

% Material boundaries (meters)
% Format: [start, end] for each material
% Each material gets a new row
material_boundaries = [0.45, 0.60;
                       0.70, 0.90;];
% generate material properties for entire domain
[eps_r, mu_r, sigma_e, sigma_m] = get_material_profile(...
    x, material_boundaries, eps_r_list, mu_r_list, sigma_e_list, sigma_m_list);

% Frequency parameters
frequency_range      = [300e6, 600e6];
frequency_resolution = 2e6;
frequencies = frequency_range(1):frequency_resolution:frequency_range(2);

% Physical constants
eps0 = 8.854187817e-12;
mu0  = 4*pi*1e-7;
c0   = 1/sqrt(eps0*mu0);

%% Boundary definitions 'o' for open or 'c' closed
left_boundary = 'o';
right_boundary = 'c';

%% Visualize regions
plot_regions(domainLength, material_boundaries, ...
             eps_r_list, mu_r_list, sigma_e_list, sigma_m_list, ...
             left_boundary, right_boundary);

%% Frequency Sweep
[Reflection_Coefficient, Transmission_Coefficient] = sweep_frequencies( ...
    frequencies, x, eps_r, mu_r, sigma_e, sigma_m, ...
    eps0, mu0, left_boundary, right_boundary);

%% Plot Results
plot_results(frequencies, Reflection_Coefficient, Transmission_Coefficient);


%% Experimentation 5-Layer Optimization 
% Goal -90dB reflection coefficient @ 500 MHz
if run_layer_optimization == true
%% 5-Layer Optimization
f0      = 500e6;
Ltot    = 0.7;     %overall lenght of structure in meters
target  = -90;      %dB                             % dB
nLayers = 5;

%lower and upper bounds 
%relative permitivity between 1.6 and 7
%electric conductivity between 0.002 and 0.0105
%thickness weights for 5 layers 0.02 to 1.00
lb = [1.6*ones(1,5), 0.002*ones(1,5),  0.02*ones(1,5)];
ub = [7.0*ones(1,5), 0.0105*ones(1,5), 1.00*ones(1,5)];

%built in MATLAB optimization routine
opts = optimoptions('lsqnonlin', 'Display', 'off', ...
    'FunctionTolerance', 1e-20, 'StepTolerance', 1e-14, ...
    'OptimalityTolerance', 1e-20, 'MaxIterations', 400);
%cost function is the layer residual for a choice of p 
%p is a 1x15 vector of free variables organized like:
%p(1:5)     relative permittivity, layers 1-5
%p(6:10)    electric conductivity (S/m), layers 1-5
%p(11:15)   thickness weights, layers 1-5
fun = @(p) layer_residual(p, x, f0, eps0, mu0, Ltot, domainLength, 'o');

best.dB = inf;
rng(1);
%try 100 times to optimize p
for trial = 1:100
    p0 = lb + rand(size(lb)).*(ub - lb);         % random start
    [p, ~] = lsqnonlin(fun, p0, lb, ub, opts);
    r  = fun(p);
    dB = 20*log10(hypot(r(1), r(2))/1e4);
    if dB < best.dB
        best.dB = dB;  best.p = p;
        fprintf('trial %3d: |Gamma| = %.1f dB\n', trial, dB);
    end
    if best.dB < target - 10, break; end         % 10 dB margin, then stop
end

%% Unpack, rebuild, and verify across frequency
p = best.p;
eps_best = p(1:5);  sig_best = p(6:10);
d_best   = Ltot*p(11:15)/sum(p(11:15));
edges    = (domainLength-Ltot)/2 + [0, cumsum(d_best(:)).'];
%% Plot optimized regions
material_boundaries_opt = [edges(1:end-1).', edges(2:end).'];   % 5x2 [start, end]

plot_regions(domainLength, material_boundaries_opt, ...
             eps_best, ones(1,5), sig_best, zeros(1,5), ...
             'o', 'o');
[eps_r, sigma_e] = build_profile(x, edges, eps_best, sig_best);
mu_r = ones(size(x));  sigma_m = zeros(size(x));
%% Rerun frequency sweep
frequencies = 300e6:2e6:600e6;
[Gam, T] = sweep_frequencies(frequencies, x, eps_r, mu_r, sigma_e, sigma_m, ...
                             eps0, mu0, 'o', 'o');
plot_results(frequencies, Gam, T);

%% Plot phase and magnitude of E_tot for this structure
[E_total, E_scattered] = calculate_fields(f0, x, eps_r, mu_r, sigma_e, sigma_m, ...
                                       eps0, mu0, 'o', 'o');
figure;
plot(x, abs(E_scattered), 'LineWidth', 1.5);
xlabel('Position (m)'); ylabel('|E_{tot}|');
title('Scattered Electric Field Magnitude at 500 MHz');
grid on;

figure;
plot(x, angle(E_scattered), 'LineWidth', 1.5);
xlabel('Position (m)'); ylabel('phase(E_{tot})');
title('Scattered Electric Field Phase at 500 MHz');
grid on;
end