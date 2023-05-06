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
z=2000; %propagation dist (m)

w_in = 1e-2;
x1c = 4e-2;
y1c = 0;
x2c = -4e-2;
y2c = 0;

[X1,Y1]=meshgrid(x1,y1);
u1 =  exp(-((X1-x1c).^2 + (Y1-y1c).^2)/(w_in^2)); %second HG00 signal
u2 =  exp(-((X1-x2c).^2 + (Y1-y2c).^2)/(w_in^2)); %second HG00 signal

uin = u1 + u2;
% u1 = u1/sqrt(int2(abs(u1).^2, sim_reg));
I1=abs(uin.^2); %src irradiance
%
figure(1)
imagesc(x1,y1,I1);
axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)');
title('z= 0 m')

uout=propTF(uin,L1,lambda,z); %propagation
%u2=propIR(u1,L1,lambda,z); %propagation
x2=x1; %obs coords
y2=y1;
I2=abs(uout.^2); %obs irrad
figure(2) %display obs irrad
imagesc(x2,y2,nthroot(I2,3));%stretch image contrast
axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)');
title(['z= ',num2str(z),' m']);
