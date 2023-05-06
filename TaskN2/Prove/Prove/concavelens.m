
L1=25e-3; %side length
M=5000; %number of samples
dx1=L1/M; %src sample interval
x1=-L1/2:dx1:L1/2-dx1; %src coords
y1=x1;
sim_reg = SimulationRegion(L1, L1, M, M);

lambda=0.5*10^-6; %wavelength
% z=2*250; %propagation dist (m)
z=500;
w_in = 0.5e-3; % Transmitter spot size
rl = 1e-3;

wr = 25e-3;
zf=0.25;
r=0.0125;

[X1,Y1]=meshgrid(x1,y1);

uin =  zeros(length(X1),length(Y1));

for i = -1:1
    for j = -1:1
        xc = rl*i;
        yc = rl*j;
        uin = uin +  exp(-((X1-xc).^2 + (Y1-yc).^2)/(w_in^2));
    end
end


u4=propIR(uin,L1,lambda,z);
% Figure input and output signals
I1=abs(u1.^2); %src irradiance
figure(1)
imagesc(x1,y1,I1);
axis square; axis xy;
% colormap('gray');  xlabel('x (m)'); ylabel('y (m)');
title('original signal');
I4=abs(u4.^2);
figure(4) %display obs irrad
imagesc(x1,y1,I4);
axis square; axis xy;
% colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title(['signal at the plane which is 10*f far from lens shows 9x magnification' ...
    '']);