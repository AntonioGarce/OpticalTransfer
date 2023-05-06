L1=0.5; %side length
M=2000; %number of samples
dx1=L1/M; %src sample interval
x1=-L1/2:dx1:L1/2-dx1; %src coords
y1=x1;

z=8000; %propagation dist (m)

%parameters for HG00 input signals
w_in = 1e-2;        %width of HG00 signal
r = 2e-2;

%generate two HG00 input signals
[X1,Y1]=meshgrid(x1,y1);

% disp_v = r;
% uin =  exp(-((X1-r).^2 + Y1.^2)/(w_in^2));
% for i = 1:8
%     disp_v = disp_v * (cos(2*pi/9)+1j*sin(2*pi/9));
%     xc = real(disp_v);
%     yc = imag(disp_v);
%     uin = uin +  exp(-((X1-xc).^2 + (Y1-yc).^2)/(w_in^2));
% end
uin = 0;
for i = -1:1
    for j = -1:1
%         if((i~=0) && (j~=0))
            xc = r*i;
            yc = r*j; 
            uin = uin +  exp(-((X1-xc).^2 + (Y1-yc).^2)/(w_in^2));
%         end
    end
end

I1=abs(uin.^2); %src irradiance

%Plot input signal 
figure(1)
imagesc(x1,y1,I1);
axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)');
title('z= 0 m');
% uout=propTF(uin,L1,lambda,z); %propagation
uout=propIR(uin,L1,lambda,z); %propagation
x2=x1; %obs coords
y2=y1;
I2=abs(uout.^2); %obs irrad

%Plot output signal
figure(2) %display obs irrad
imagesc(x2,y2,nthroot(I2,3));%stretch image contrast
axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)');
title(['z= ',num2str(z),' m']);