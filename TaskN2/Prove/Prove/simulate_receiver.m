tartup;

%% Simulation parameters

% Define the simulation window
sim_reg = SimulationRegion(10e-2, 10e-2, 512, 512);

% Define channel properties
L = 4;
n_screen=1; % Number of phase screen

% Input beam properties 
w_in = 7e-3/2; % Transmitter spot size
U_in= exp(-(sim_reg.X.^2 + sim_reg.Y.^2)/(w_in^2)); % normalized gaussian beam used as input of the turbulent channel
U_in = U_in/sqrt(int2(abs(U_in).^2, sim_reg));

% Turbulence parameters (Modified Von Karman)
Cn2=1e-12; % [m^(-2/3)]
D_0=1000; % Outer scale   [m]
d_0=1e-6; %  Inner scale  [m]
corr_coeff = 0.9; % temporal correlation coefficient

% Simulation parameters
n_iter = 200; % Iteration number

% Receiver properties
S = 2e-2;
receiver = Receiver2(sim_reg, [S/2, 0]);
receiver.w0 = 3e-3/2;
receiver.is_symmetric = true;

%% Simulate turbulent propagation
U_out = zeros(sim_reg.N_x, sim_reg.N_y, n_iter);
phase_no_turb = angle(free_prop(U_in, sim_reg, L));
U_rec = receiver.U.*exp(1i*phase_no_turb);
phase_screen = zeros(sim_reg.N_x, sim_reg.N_y, n_iter);
c = zeros(1,n_iter);
phase_screen_old=zeros(sim_reg.N_x, sim_reg.N_y);

for i = 1:n_iter
    % propagation through a turbulent medium
    [U_out(:, :, i), phase_screen_old] = turbulent_prop(U_rec, sim_reg, L, Cn2, D_0, d_0, n_screen, phase_screen_old, corr_coeff);
    phase_screen(:, :, i) = phase_screen_old;
    
    % overlap between the incoming and the receiver's fields
    c(i) = overlap(receiver.U ,U_out(:, :, i), sim_reg);
end
%% Plot the receiver geometry
[fig, ax] = plot_receiver(receiver, sim_reg, S, receiver.w0, r0(Cn2, L, k));

%% Post-process
M(n_iter) = struct('cdata',[],'colormap',[]);

fig = figure();
fig.Visible = "on";
tiledlayout(1, 1, "TileSpacing", "compact");
ax = nexttile;

for i = 1:n_iter
    img = imagesc(ax, sim_reg.x*1e3, sim_reg.y*1e3, abs(conj(U_rec).*U_out(:, :, i)).^2);
    title('Captured intensity at Receiver');
    xlabel('x [mm]');
    ylabel('y [mm]');
    axis(ax, "square");
    colormap(ax, "hot")
    drawnow;
    M(i) = getframe;
end   
%% Plot normalized intensity
eta = abs(c).^2;
eta_norm = eta - mean(eta);

fig = figure();
hold on;
plot(eta, "DisplayName", "no normalization");
plot(eta_norm, "DisplayName", "normalized");
xlabel("iteration number");
ylabel("coupling efficiency");
grid("on");
box("on");
legend();