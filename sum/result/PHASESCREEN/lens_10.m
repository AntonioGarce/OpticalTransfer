% sqr_beam propagation example
%
clear;
close("all");
clc;
startup;


L1=25e-3; %side length
M=5000; %number of samples
dx1=L1/M; %src sample interval
x1=-L1/2:dx1:L1/2-dx1; %src coords
y1=x1;
sim_reg = SimulationRegion(L1, L1, M, M);

lambda=0.5*10^-6; %wavelength
k=2*pi/lambda; %wavenumber
w=0.01; %source half width (m)
z=2*250; %propagation dist (m)

w_in = 2.5e-3; % Transmitter spot size

wr = 10e-3;
window=circ(x1,y1,wr);

L = 100;
Cn2=1e-10;
D_0=1000; % Outer scale   [m]
d_0=1e-6; %  Inner scale  [m]
corr_coeff = 0.1;
% Simulation parameters
n_screen=1; % Number of phase screen
n_iter = 10; % Iteration number
% Generate a normalized gaussian beam as input
u1= exp(-(sim_reg.X.^2 + sim_reg.Y.^2)/(w_in^2)); % normalized gaussian beam used as input of the free channel
u1 = u1/sqrt(int2(abs(u1).^2, sim_reg));
phase_screen_old=zeros(sim_reg.N_x, sim_reg.N_y);
[u1(:, :), phase_screen_old] = turbulent_prop(u1, sim_reg, L, Cn2, D_0, d_0, n_screen, phase_screen_old, corr_coeff);
% phase_screen(:, :) = phase_screen_old;
zf=0.25;
r=0.0125;
u1=window.*u1;
P=circ(x1,y1,r);
u1=P.*u1;

I1=abs(u1.^2); %src irradiance

figure(1)
imagesc(x1,y1,I1);
axis square; axis xy;
% colormap('gray');  xlabel('x (m)'); ylabel('y (m)');
title(['original signal']);

u1=focus(u1,L1,lambda,zf);
u5=propIR(u1,L1,lambda,4*zf);
u2=propTF(u1,L1,lambda,2*zf);
u2=P.*u2;
u3=focus(u2,L1,lambda,zf);
u4=propTF(u3,L1,lambda,3*zf);


x2=x1; %obs coords
y2=y1;


I2=abs(u2.^2); %obs irrad
figure(2) %display obs irrad
imagesc(x2,y2,I2);
axis square; axis xy;
% colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title(['signal at the second lens plane']);

I4=abs(u4.^2);
figure(3) %display obs irrad
imagesc(x1,y1,I4);
axis square; axis xy;
% colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title(['final signal after the second lens']);

I5=abs(u5.^2);
figure(4) %display obs irrad
imagesc(x1,y1,I5);
axis square; axis xy;
% colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title(['final signal after the second lens']);


