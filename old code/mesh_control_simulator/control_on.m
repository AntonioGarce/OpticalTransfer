function [output, step_colormap, heater1, heater2] = control_on(in,N,kp0,ag,freq,fC,fC_idx,lambdaC,c0,neff0,ng,Lmzi,delta_L,H1,H2,L)

% define the output matrices to be filled
output = zeros(N+1,length(freq));
step_colormap = zeros(length(H1), length(H2), N);
heater1 = zeros(1,N);
heater2 = zeros(1,N);

t11 = zeros(1,length(freq));
t12 = zeros(1,length(freq));
t21 = zeros(1,length(freq));
t22 = zeros(1,length(freq));

T11 = zeros(N,length(freq));
T12 = zeros(N,length(freq));
T21 = zeros(N,length(freq));
T22 = zeros(N,length(freq));

input = zeros(1,2);
ol_compensation = zeros(1,length(freq));  % optical length compensation


input(1,1) = in(1);  % take first two inputs
input(1,2) = in(2);
in2 = in(2);

stepcross = input(1,1);


for conta=1:N
    
    disp(['Starting control on step ', num2str(conta), ' . . .']);
    
    [colormap, h1_sel, h2_sel]=heaters_selection(input,kp0,ag,fC,neff0,Lmzi,delta_L,H1,H2,L,conta); % fa il controllo sul primo step
    % save results
    step_colormap(:,:,conta) = colormap;
    heater1(conta) = h1_sel*pi;
    heater2(conta) = h2_sel*pi;
    
    for it=1:length(freq)
        
        f=freq(it);
        
        neff=neff0+((ng-neff0)/fC)*(f-fC);    % neff dispersion
        
        slope=0.015/(40e-9);
        kp=kp0-((lambdaC^2)*(f-fC)*slope)/c0; % k dispersion
        
        [mz]=MZI(kp,ag,f,neff,Lmzi,delta_L,heater1(conta),heater2(conta),L,conta);
        
        t11(1,it) = mz(1,1);
        t12(1,it) = mz(1,2);
        t21(1,it) = mz(2,1);
        t22(1,it) = mz(2,2);
        
        ol_compensation(1,it) = exp((ag -1i*2*pi*f*neff/c0) * Lmzi * conta);
        
    end
    
    T11(conta,:) = t11;
    T12(conta,:) = t12;
    T21(conta,:) = t21;
    T22(conta,:) = t22;
    
    % lower output
    output(conta,:) = abs( T11(conta,:).*stepcross+T12(conta,:).*in2 ).^2;
    % upper output
    stepcross = T21(conta,:).*stepcross+T22(conta,:).*in2;
    output(conta+1,:) = abs(stepcross).^2;
    
    
    % Check control
    bar_minimized = output(conta,:);
    cross_maximized = output(conta+1,:);
    disp(['Bar minimized value: ', num2str(bar_minimized(1,fC_idx))]);
    disp(['Cross maximized value: ', num2str(cross_maximized(1,fC_idx))]);
    disp(' ');
    
    
    
    % set input for next step control
    input(1,1) = stepcross(1,fC_idx);
    in2 = in(conta+2) .* ol_compensation;
    input(1,2) = in2(1,fC_idx);
    
    
end




end

