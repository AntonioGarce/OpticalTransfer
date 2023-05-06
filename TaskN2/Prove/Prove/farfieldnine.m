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
z0=0.0625;
zf=0.25;
r=2e-3;


k=2*pi/lambda; 

disp_v = r;

X1=sim_reg.X;
Y1=sim_reg.Y;

uin =  exp(-((X1-r).^2 + Y1.^2)*(1/(w_in^2)+1i*k/(2*z0)));

for i = 1:8
    disp_v = disp_v * (cos(2*pi/9)+1j*sin(2*pi/9));
    xc = real(disp_v);
    yc = imag(disp_v);
    uin = uin +  exp(-((X1-xc).^2 + (Y1-yc).^2)*(1/(w_in^2)+1i*k/(2*z0)));
end


u1 = uin/sqrt(int2(abs(uin).^2, sim_reg));

% P=circle_window(x1,y1,wr);
% u1=P.*u1;


u2=propTF(u1,L1,lambda,zf-z0);
u2_dot=focus(u2,L1,lambda,zf);

u3=propTF(u2_dot,L1,lambda,100);
u3_dot=focus(u3,L1,lambda,zf);

uout=propTF(u3_dot,L1,lambda,zf+z0);

I1=abs(uin.^2); %src irradiance

figure(1);
imagesc(x1,y1,I1);
axis square; axis xy;
title('original signal:hg00');

I2=abs(u2.^2);
figure(2) %display obs irrad
imagesc(x1,y1,I2);
axis square; axis xy;
title('farfield signal');

I3=abs(uout.^2);
figure(3) %display obs irrad
imagesc(x1,y1,I3);
axis square; axis xy;
title('output signal');

