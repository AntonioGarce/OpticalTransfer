%Set up the workspace
clc;clear
startup;

%% Simulation parameters

% Define the simulation window
sim_reg = SimulationRegion(5e-2, 5e-2, 512, 512);
% Define channel properties
L = 7;
Cn2=1e-10;
D_0=1000; % Outer scale   [m]
d_0=1e-6; %  Inner scale  [m]
corr_coeff = 0.0001;
% Simulation parameters
n_screen=1; % Number of phase screen
n_iter = 1; % Iteration number

% Transmitter properties 
w_in = 1e-2; % Transmitter spot size

% Generate a normalized gaussian beam as input

U_in= exp(-(sim_reg.X.^2 + sim_reg.Y.^2)/(w_in^2)); % normalized gaussian beam used as input of the free channel
U_in = U_in/sqrt(int2(abs(U_in).^2, sim_reg));
xc = 0;
yc = 0; % initial coordinates of aperture
xcount = 0; %counter to cicle through positions of the aperture
d = 50e-6;
c_free = zeros(9,1);% Preallocate memory

c_turbulent = zeros(n_iter,9);

% % disp = [5e-6, 5e-6];
% % receiver = Receiver9(sim_reg,disp);
% % U_out = zeros(sim_reg.N_x, sim_reg.N_y, n_iter);
% % phase_no_turb = angle(free_prop(U_in, sim_reg, L));
% % U_rec = receiver.U.*exp(1i*phase_no_turb);
% % phase_screen = zeros(sim_reg.N_x, sim_reg.N_y, n_iter);
% % phase_screen_old=zeros(sim_reg.N_x, sim_reg.N_y);
% % 
% % [U_out(:, :), phase_screen_old] = turbulent_prop(U_rec, sim_reg, L, Cn2, D_0, d_0, n_screen, phase_screen_old, corr_coeff);
% % phase_screen(:, :) = phase_screen_old;
% % % overlap between the incoming and the receiver's fields
% % c_turbulent = overlap(receiver.U ,U_out(:, :), sim_reg);


for k = 1:9
% Receiver properties
    S = 2e-2;
    
    receiver = Receiver1(sim_reg);
    receiver.w0 = 3e-3/2;
    receiver.center = [xc, yc]; 
    

    U_out = free_prop(U_in, sim_reg, L);
%     U = imagesc(receiver.U);
    c = overlap(receiver.U ,U_out(:, :), sim_reg);
    c_free(k) = c;

    if xcount<2
        xc=xc+d;
        xcount=xcount+1;

    else
        xc=0;
        xcount=0;
        yc=yc+d;
    end

% fprintf("complex overlap: %f \n", c);
% fprintf("power overlap [dB]: %f \n", 10*log10(abs(c)^2));

%% Simulate turbulent propagation
    U_out = zeros(sim_reg.N_x, sim_reg.N_y, n_iter);
    phase_no_turb = angle(free_prop(U_in, sim_reg, L));
    U_rec = receiver.U.*exp(1i*phase_no_turb);
    phase_screen = zeros(sim_reg.N_x, sim_reg.N_y, n_iter);
    phase_screen_old=zeros(sim_reg.N_x, sim_reg.N_y);

  
    % propagation through a turbulent medium
        [U_out(:, :), phase_screen_old] = turbulent_prop(U_rec, sim_reg, L, Cn2, D_0, d_0, n_screen, phase_screen_old, corr_coeff);
        phase_screen(:, :) = phase_screen_old;
    % overlap between the incoming and the receiver's fields
        c_turbulent(k) = overlap(receiver.U ,U_out(:, :), sim_reg);
    
end
I1=abs(U_in.^2); %src irradiance
figure(1)
imagesc(sim_reg.x,sim_reg.y,I1);
axis square; axis xy;
%     colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title('original signal')
fnl=num2str(k)+"in.fig";
saveas(gcf,fnl);

I2=abs(U_out.^2); %src irradiance
figure(2)
imagesc(sim_reg.x,sim_reg.y,I2);
axis square; axis xy;
%     colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title('signal after phase screen')
fnl="out.fig";
saveas(gcf,fnl);

% % U_last = receiver.U.*U_out;
% % Il=abs(U_last.^2); %src irradiance
% % figure(3);
% % imagesc(sim_reg.x,sim_reg.y,Il);
% % saveas(gcf,"received.fig");
c_free=conj(c_free);
c_free=c_free/norm(c_free);
c_turbulent=c_turbulent/norm(c_turbulent);
orthogonality=-dot(c_free,c_turbulent);
