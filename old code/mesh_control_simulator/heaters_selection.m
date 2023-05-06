% This code is to check the control for one MZI step. 
% It creates for each step the colormap for the port that must be minimized. 
% Then, we search for the min in the colormap and returns the step phases (H1 and H2), 
% which are associated to the heaters.


function [colormap, h1_selected, h2_selected]=heaters_selection(in,kp0,ag,fC,neff0,Lmzi,delta_L,H1,H2,L,conta)


% create the colormap

colormap = zeros(length(H1), length(H2), 1);

for i = 1 : length(H1)
    
    h1 = H1(i);
%     disp(['Pre-step heater phase: ', num2str(h1/pi), 'pi']);
    
    for j = 1 : length(H2)
        
        h2 = H2(j);
        
        [mz] = MZI(kp0,ag,fC,neff0,Lmzi,delta_L,h2,h1,L,conta);   % only evaluated at fC
        
        t11 = mz(1,1);
        t12 = mz(1,2);
        t21 = mz(2,1);
        t22 = mz(2,2);
        
        bar = abs(  t11*in(1)+t12*in(2)  ).^2;
        cross = abs(  t21*in(1)+t22*in(2)  ).^2;

        colormap(i,j,1) = bar;
        
    end
end

% find the minimum and isolate h1 and h2 phases

tmp = colormap;
minimum = min(min(tmp(:)));
[x,y,~]=find(tmp==minimum, 1, 'last');

h1_selected = H1(y)/pi;
h2_selected = H2(x)/pi;

disp(' ');
disp(['Selected H1: ', num2str(h1_selected), 'pi']);
disp(['Selected H2: ', num2str(h2_selected), 'pi']);




end