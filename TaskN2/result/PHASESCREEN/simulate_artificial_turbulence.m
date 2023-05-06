% Set up the workspace
startup;

%% Simulation parameters

% Define the simulation window
sim_reg = SimulationRegion(5e-2, 5e-2, 512, 512);

% Turbulence parameters (Modified Von Karman)
Cn2=1e-10; % [m^(-2/3)]
D0=1000; % Outer scale   [m]
d0=1e-6; %  Inner scale  [m]
corr_coeff = 0.9; % temporal correlation coefficient

% Define channel properties
L1 = 4 - 0.125; %
L2 = 0.25; 
L3 = 4 -0.125;
L = L1 + L2 + L3;

% Simulation parameters
n_screen=1; % Number of phase screen
n_iter = 200; % Iteration number

% Transmitter properties 
w_in = 40e-3/2; % Transmitter spot size

% Receiver properties
S = 2e-2;
receiver = Receiver1(sim_reg);
receiver.w0 = 3e-3/2;

%% Simulate free propagation
% Generate a normalized gaussian beam as input
U_in= exp(-(sim_reg.X.^2 + sim_reg.Y.^2)/(w_in^2)); % normalized gaussian beam used as input of the turbulent channel
U_in = U_in/sqrt(int2(abs(U_in).^2, sim_reg));

% Preallocate memory
% phase_no_turb = angle(free_prop(U_in, sim_reg, L));
% U_rec = receiver.U.*exp(1i*phase_no_turb);

U_out = free_prop(U_in, sim_reg, L);
c = overlap(receiver.U ,U_out(:, :), sim_reg);

fprintf("complex overlap: %f \n", c);
fprintf("power overlap [dB]: %f \n", 10*log10(abs(c)^2));

%% Simulate turbulent propagation
% Generate a normalized gaussian beam as input
U_in= exp(-(sim_reg.X.^2 + sim_reg.Y.^2)/(w_in^2)); % normalized gaussian beam used as input of the turbulent channel
U_in = U_in/sqrt(int2(abs(U_in).^2, sim_reg));

% Preallocate memory
U_out = zeros(sim_reg.N_x, sim_reg.N_y, n_iter);
phase_no_turb = angle(free_prop(U_in, sim_reg, L));
U_rec = receiver.U.*exp(1i*phase_no_turb);
phase_screen = zeros(sim_reg.N_x, sim_reg.N_y, n_screen, n_iter);
c = zeros(1,n_iter);
phase_screen_old=zeros(sim_reg.N_x, sim_reg.N_y);

tic;

U_in = free_prop(U_in, sim_reg, L1);
for i = 1:n_iter
    % propagation through a concentrated turbulent medium
    [U_out(:, :, i), phase_screen_old] = turbulent_prop(U_in, sim_reg, L2, Cn2, D0, d0, n_screen, phase_screen_old, corr_coeff);
    phase_screen(:, :, :, i) = phase_screen_old;
    
    U_out(:, :, i) = free_prop(U_out(:, :, i), sim_reg, L3);
    % overlap between the incoming and the receiver's fields
    c(i) = overlap(receiver.U ,U_out(:, :, i), sim_reg);
end

computation_time = toc;

fprintf("computation time [s]: %f \n", computation_time)

%% Phase screen number sweep
n_screen_sweep = 1:10;

% Generate a normalized gaussian beam as input
U_in= exp(-(sim_reg.X.^2 + sim_reg.Y.^2)/(w_in^2)); % normalized gaussian beam used as input of the turbulent channel
U_in = U_in/sqrt(int2(abs(U_in).^2, sim_reg));

U_out = cell(size(n_screen_sweep, 1), size(n_screen_sweep, 2));
phase_screen = cell(size(n_screen_sweep, 1), size(n_screen_sweep, 2));

tic;

for n = 1:length(U_out)
    U_out{n} = zeros(sim_reg.N_x, sim_reg.N_y, n_iter);
    phase_screen{n} = zeros(sim_reg.N_x, sim_reg.N_y, n_screen_sweep(n), n_iter);
    phase_screen_old=zeros(sim_reg.N_x, sim_reg.N_y);
    
    U_L1 = free_prop(U_in, sim_reg, L1);
    for i = 1:n_iter
        % propagation through a concentrated turbulent medium
        [U_out{n}(:, :, i), phase_screen_old] = turbulent_prop(U_L1, sim_reg, L2, Cn2, D0, d0, n_screen_sweep(n), phase_screen_old, corr_coeff);
        phase_screen{n}(:, :, :, i) = phase_screen_old;
        U_out{n}(:, :, i) = free_prop(U_out{n}(:, :, i), sim_reg, L3);
    end
end
clear("U_L1");

computation_time = toc;

fprintf("computation time [s]: %f \n", computation_time)

%% Cn2 sweep
Cn2_sweep = 10.^(-15:-9);

% Generate a normalized gaussian beam as input
U_in= exp(-(sim_reg.X.^2 + sim_reg.Y.^2)/(w_in^2)); % normalized gaussian beam used as input of the turbulent channel
U_in = U_in/sqrt(int2(abs(U_in).^2, sim_reg));

U_out = cell(size(Cn2_sweep, 1), size(Cn2_sweep, 2));
phase_screen = cell(size(Cn2_sweep, 1), size(Cn2_sweep, 2));

tic;

for n = 1:length(Cn2_sweep)
    U_out{n} = zeros(sim_reg.N_x, sim_reg.N_y, n_iter);
    phase_screen{n} = zeros(sim_reg.N_x, sim_reg.N_y, n_iter);
    phase_screen_old=zeros(sim_reg.N_x, sim_reg.N_y);
    
    U_L1 = free_prop(U_in, sim_reg, L1);
    for i = 1:n_iter
        % propagation through a concentrated turbulent medium
        [U_out{n}(:, :, i), phase_screen_old] = turbulent_prop(U_L1, sim_reg, L2, Cn2_sweep(n), D0, d0, n_screen, phase_screen_old, corr_coeff);
        phase_screen{n}(:, :, i) = phase_screen_old;
        U_out{n}(:, :, i) = free_prop(U_out{n}(:, :, i), sim_reg, L3);
    end
end
clear("U_L1");

computation_time = toc;

fprintf("computation time [s]: %f \n", computation_time)

%% Compute effective Cn2 from scintillation index
% Define the integration window for the scintillation index (TO BE
% COMPLETED)
win_width = 10; % window length in px
win_heigth = 10; % window heigth in px

sigma_I = zeros(size(Cn2_sweep, 1), size(Cn2_sweep, 2));
for n = 1:size(Cn2_sweep, 2)
    tmp = scint_idx(abs(U_out{n}).^2);
    sigma_I(n) = tmp(sim_reg.N_x/2, sim_reg.N_y/2);
end
clear tmp;
Cn2_eff = Cn(sigma_I, k, L);

save( ...
    "data\Cn2_effective_10m_C.mat", ...
    "Cn2_sweep", ...
    "Cn2_eff", ...
    "sigma_I", ...
    "sim_reg",...
    "Cn2", ...
    "D0", ...
    "d0", ...
    "corr_coeff", ...
    "L1", ...
    "L2", ...
    "L3", ...
    "L", ...
    "n_screen", ...
    "n_iter", ...
    "w_in" ...
);

%% Plot the relationship between Cn2 and the Cn2 effective
fig = figure();
tiledlayout(1, 1);

ax = nexttile();
loglog(Cn2_sweep, Cn2_eff);

title("Effective $C_n^2$ vs real $C_n^2$");
xlabel("$C_n^2$");
ylabel("$C_{n, eff}^2$");
axis("tight");
grid("on");
box("on");

%% Plot the comparison between r0 real and corrected
ax = nexttile();
hold on;
plot(Cn2_eff, r0(Cn2_eff, k, L), "DisplayName", "effective $r_0$");
plot(Cn2_eff, r0(Cn2_sweep, k, L2), "DisplayName", "corrected $r_0$");
hold off;

set(ax, "XScale", "log");
set(ax, "YScale", "log");
xlabel("$C_{n, eff}^2$");
ylabel("Fried Parameter");
axis("tight");
legend();
grid("on");
box("on");

%% Plot a random realization of the output beam
close all
plot_angle(U_out(:, :, randi(n_iter)), sim_reg)
[fig, ax, foo] = plot_int(U_out(:, :, randi(n_iter)), sim_reg, "Normalize", true);
title(ax, "Simulated intensity at the receiver")

%% Plot scintillation index
sigma_I = scint_idx(abs(U_out).^2);
sigma_I(sigma_I > 1e-2) = nan;
plot_int(sqrt(sigma_I), sim_reg, "normalize", false);

fprintf("%f \n", sigma_I(sim_reg.N_x/2, sim_reg.N_y/2))
fprintf("%f \n", sim_reg.X(sim_reg.N_x/2, sim_reg.N_y/2))
fprintf("%f \n", sim_reg.Y(sim_reg.N_x/2, sim_reg.N_y/2))
%% Estimate the correlation radius of the output beam
std_phase_data = std_phase(U_out, sim_reg, [sim_reg.N_x/2, sim_reg.N_y/2]);

r0 = regionprops(std_phase_data < 1, "MajorAxisLength");

corr_region = std_phase_data;
corr_region(corr_region > 1) = nan; 
corr_region = exp(1i*corr_region);

[fig, ax] = plot_int(U_out(:, :, end), sim_reg)

[fig, ax, cbar] = plot_angle(corr_region, sim_reg);
title(ax, "Correlated region");
cbar.Limits = [0, 1];
cbar.Ticks = 0:0.1:1;
cbar.TickLabelsMode = "auto";
cbar.Label.String = "$\sigma(\Delta\phi)$";
cbar.Label.Interpreter = "latex";

%% Evaluate DOC
DOC_centered = zeros(sim_reg.N_x, sim_reg.N_y);
for i = 1:sim_reg.N_x
    for j = 1:sim_reg.N_y
        DOC_centered(i, j) = DOC(U_out, [i, j], [0, 0]);
    end
end

imagesc(sim_reg.X, sDOC_centered)