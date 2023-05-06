L1=25e-3; %side length
M=5000; %number of samples
dx1=L1/M; %src sample interval
x1=-L1/2:dx1:L1/2-dx1; %src coords
y1=x1;
sim_reg = SimulationRegion(L1, L1, M, M);

lambda=0.5*10^-6; %wavelength
% z=2*250; %propagation dist (m)

w_in = 0.5e-3; % Transmitter spot size

wr = 25e-3;
z0=0.025;
zf=0.25;
r=0.0125;

x1c = 1e-3;         %x-axis position of center of first HG00 signal
y1c = 0;            %y-axis position of center of first HG00 signal
x2c = -1e-3;        %x-axis position of center of second HG00 signal
y2c = 0;            %y-axis position of center of second HG00 signal
k=2*pi/lambda; 
uin1= exp(-((sim_reg.X-x1c).^2 + (sim_reg.Y-y1c).^2)*(1/(w_in^2)+1i*k/(2*z0)));
uin2=exp(-((sim_reg.X-x2c).^2 + (sim_reg.Y-y2c).^2)*(1/(w_in^2)+1i*k/(2*z0))); % normalized gaussian beam used as input of the free channel

u1=uin1+uin2;
u1 = u1/sqrt(int2(abs(u1).^2, sim_reg));

P=circle_window(x1,y1,r);
u1=P.*u1;

% uin1_div=uin1.*exp(1i*k/(2*z0)*((sim_reg.X-x1c).^2+(sim_reg.Y-y1c).^2)); %diverging hg00
% uin2_div=uin2.*exp(1i*k/(2*z0)*((sim_reg.X-x2c).^2+(sim_reg.Y-y2c).^2)); %diverging hg00
% u1_div=uin1_div+uin2_div;

u2=propTF(u1 ,L1,lambda,zf-z0);
u2_dot=focus(u2,L1,lambda,zf);

u3=propTF(u2_dot,L1,lambda,zf);
u3_dot=focus(u3,L1,lambda,zf);

uout=propTF(u3_dot,L1,lambda,zf-z0);

I1=abs(u1.^2); %src irradiance

figure(1);
imagesc(x1,y1,I1);
axis square; axis xy;
title('original signal:hg00');

I2=abs(u2.^2);
figure(2) %display obs irrad
imagesc(x1,y1,I2);
axis square; axis xy;
title('input signal-diverging hg00');

I3=abs(uout.^2);
figure(3) %display obs irrad
imagesc(x1,y1,I3);
axis square; axis xy;
title('input signal-diverging hg00');

