%Set up the workspace
clc;clear
startup;

L1=0.5; %side length
M=2000; %number of samples
dx1=L1/M; %src sample interval
x1=-L1/2:dx1:L1/2-dx1; %src coords
y1=x1;
lambda=0.5*10^-6; %wavelength
k=2*pi/lambda; %wavenumber
w=0.051; %source half width (m)
z=200; %propagation dist (m)

w_in = 1e-2;

[X1,Y1]=meshgrid(x1,y1);
u1= exp(-(X1.^2 + Y1.^2)/(w_in^2)); % normalized gaussian beam used as input of the free channel
% u1 = u1/sqrt(int2(abs(u1).^2, sim_reg));
I1=abs(u1.^2); %src irradiance
%
figure(1)
imagesc(x1,y1,I1);
axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)');
title('z= 0 m')
deg=pi/180;
alpha=5.0e-5; %rad
theta=45*deg;
[u1]=tilt(u1,L1,lambda,alpha,theta);
u2=propTF(u1,L1,lambda,z); %propagation
%u2=propIR(u1,L1,lambda,z); %propagation
x2=x1; %obs coords
y2=y1;
I2=abs(u2.^2); %obs irrad
figure(2) %display obs irrad
imagesc(x2,y2,nthroot(I2,3));%stretch image contrast
axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)');
title(['z= ',num2str(z),' m']);
