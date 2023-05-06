% Code to simulate a mesh made by 8 balanced mach-zehnders in cascade.
% Mesh is 9x1. Heaters are located above each MZI upper arm and before each
% MZI in the lower arm.
%
%  [9]                               x                           [1]->OUT
%       .                                                        .
%       .                                                      .
%       .                                                  (8 steps) 
%       .                                                   .
%       ________________________________/--heater--\_____ .
%       __________/--heater--\__heater__   (step2)  _____
%                    (step1)            \----------/
%  IN-> --heater--\__________/---------------------------
%
%
%
% 


clear all
close all
clc


% simulation parameters

c0=299792458;
ng = 4.2;
neff0 = 2.4;

lambda2=1570e-9;     % simulation band
lambda1=1530e-9;

deltaWL=lambda2-lambda1;
lambdaC=lambda1+deltaWL/2;
fC=c0/lambdaC;

f1=c0/lambda1;
f2=c0/lambda2;



N = 8;      % number of MZIs in the mesh

cases = 1;
result = zeros(1,cases);  % vector containing the value of the maximized power at fC for each case study
mse = zeros(1,cases);



% MZI/mesh design

kp = 0.45;              % power coupling ratio
Lmzi = 80e-6;
L = 0;                  % physical length of the lower arm line with the heater
delta_L = 100e-6*0;
FSR = c0/(ng*delta_L);

dbovercm = -1;         % Att per cm
dboverm=dbovercm.*100; % ATT per meter
ag=dboverm./8.686;     % converte dB/m in Neper


% generate 1000 different phase values from pi to 3pi (to be associated to H1 and H2)

step = 2.5*pi/999;

H1 = pi : step : 2.5*pi;
H2 = H1;


% axes settings for plot

Resolution = 250e6;
freq = f2 : Resolution : f1;
% freq = fC;

[~, fC_idx] = min(abs(freq-fC));
% fC_idx = 1;            % use only when control is for one single frequency!

lambda_axis = c0./(freq);

fth=(freq-fC)/1e12;               % frequency ax
lax = (lambda_axis-lambdaC)*1e9;  % wavelength ax



for caseno=1:cases  % case number
    
    disp(['CASE: ', num2str(caseno), ' STARTED !']);
    
    
    % define random input vector
    
%     in = [0.5+1i*2  0.2+1i*0.5  0.3+1i*0.4  0.8+1i*0.1];
    in = rand(1,N+1)+1i*rand(1,N+1);
    in = in/sqrt(sum(abs(in).^2));
    in = [in 0];
    

    
    % check the mesh control
    
    [output,step_colormap,heater1,heater2] = control_on(in,N,kp,ag,freq,fC,fC_idx,lambdaC,c0,neff0,ng,Lmzi,delta_L,H1,H2,L);
    
    output = 10*log10(output);
    
    heater1 = heater1./pi;
    heater2 = heater2./pi;
    
    
    % to visualize in wavelength
    len_o = length(output(:,1));
    for i=1:len_o
        
        output(i,:)=flip(output(i,:));
        
    end
    
    
    
    %%%%%%% FOR SINGLE FREQ fC
    
%     result(1,caseno) = output(N+1, fC_idx);
    
    
    
    % mse calculation
    
%     tot_power = sum(abs(in).^2);
%     res = result(caseno);
%     res = 10^(res/10);     % convert in linear 
%     
%     mse(1,caseno) = MSE(caseno, res, tot_power);
    
end










%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot colormap

for j=1:N
    figure();
    imagesc(H1/pi,H2/pi,step_colormap(:,:,j)); colorbar; hold on;
    plot(heater1(j), heater2(j), 'x', 'Color', 'r', 'LineWidth', 2, 'MarkerSize', 10);
    label(j) = "Colormap for step "+j+"";
    title(label(j));
    set(gca, 'Ydir', 'normal');
    xlabel('H1 phase');
    ylabel('H2 phase');
    set(gca,'fontsize',16);
end


% plot mesh outputs

figure();
for i=1:N+1

    p(i) = plot(lax, output(i,:), 'LineWidth', 1);
    if N==1
        label(1) = "bar";
        label(2) = "cross";
    else
        label(i) = "output"+i+"";
    end
    hold on;
    grid on;

end
title('Mesh control');
xlabel('Wavelength [nm]');
ylabel('Mesh transfer function [dB]');
legend(label);
set(gca,'Color',[0.8 0.8 0.8]);
set(gca,'fontsize',16);
axis([min(lax) max(lax) min(min(output))-10 max(max(output))+10]);







%%%%%%%% PLOTS for single frequency fC %%%%%%%%%

% Plot maximized power at the last mesh output (only value at fC)

% figure();
% for i=1:cases
%     scatter(i,result(1,i),'.','b'); hold on;
% end
% title('Maximized mesh output at fC');
% xlabel('Case number');
% ylabel('Output power [dB]');
% set(gca,'Color',[0.8 0.8 0.8]);
% set(gca,'fontsize',16);
% axis([ 0 cases+1 min(min(result))-2 max(max(result))+1]);


% plot mse

% figure();
% histogram(mse);
% title('MSE histogram');
% xlabel('MSE');
% ylabel('Cases');
% set(gca,'fontsize',16);












