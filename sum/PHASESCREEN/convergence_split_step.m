%%
startup

%% Simulation parameters
% Define the simulation window
sim_reg = SimulationRegion(5e-2, 5e-2, 512, 512);

% Define channel properties
L = 8; % Channel length [m]

% Input beam properties 
w_in = 7e-3/2; % Transmitter spot size
U_in= exp(-(sim_reg.X.^2 + sim_reg.Y.^2)/(w_in^2)); % normalized gaussian beam used as input of the turbulent channel
U_in = U_in/sqrt(int2(abs(U_in).^2, sim_reg));

% Turbulence parameters (Modified Von Karman)
Cn2=1e-11; % [m^(-2/3)]
D_0=1000; % Outer scale   [m]
d_0=1e-6; %  Inner scale  [m]
corr_coeff = 0.9; % temporal correlation coefficient

% Simulation parameters
n_iter = 200; % Iteration number

n_screen_sweep = 1:5;

%%

tic;

% Preallocate memory
N_x = sim_reg.N_x;
N_y = sim_reg.N_y;
U_out = cell(1, n_screen);
phase_screen = cell(1, n_screen);
parfor n_screen = n_screen_sweep
    % Preallocate memory
    dz = L/n_screen;
    U_out{n_screen} = zeros(N_x, N_y, n_iter);
    phase_screen_old=zeros(N_x, N_y, n_screen);
    phase_screen{n_screen} = zeros(N_x, N_y, n_iter, n_screen);
    disp(n_screen);
    for i = 1:n_iter
        % propagation through a turbulent medium
        [U_out{n_screen}(:, :, i), phase_screen_old] = turbulent_prop(U_in, sim_reg, L, dz, D_0, d_0, Cn2, n_screen, phase_screen_old, corr_coeff);
        phase_screen{n_screen}(:, :, i, :) = phase_screen_old;
    end
end

computation_time = toc;

fprintf("computation time [s]: %f", computation_time);

%% Analyse scintillation index
fig = cell(1, size(n_screen_sweep, 2));
ax = cell(1, size(n_screen_sweep, 2));
for n_screen = n_screen_sweep
    fig{n_screen} = figure("Name", sprintf("n_screen: %d", n_screen), "NumberTitle", "off");
    ax{n_screen} = subplot(1, 1, 1);
    imagesc(sim_reg.x*1e3, sim_reg.y*1e3, scint_idx(abs(U_out{n_screen}).^2));
    xlabel("x [mm]");
    ylabel("y [mm]");
    axis("square");
    colormap(ax{n_screen}, "hot");
    colorbar(ax{n_screen});
end

figure();
ax = subplot(1, 1, 1);
sigma_I = zeros(1, size(n_screen_sweep, 2));
for n_screen = n_screen_sweep
    tmp = scint_idx(abs(U_out{n_screen}).^2);
    sigma_I(n_screen) = tmp(sim_reg.N_x/2, sim_reg.N_y/2);
end
plot(n_screen_sweep, sigma_I)
xlabel("screen number");
ylabel("scintillation index");
grid("on");
box("on");
%%
n_screen = 3;
n_iter = 90;

fig = figure(1);
img = imagesc(sim_reg.x*1e3, sim_reg.y*1e3, abs(U_out{n_screen}(:, :, n_iter)).^2);
xlabel("x [mm]");
ylabel("y [mm]");
axis("square");

