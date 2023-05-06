function[phase_screen]=phase_screen_MVKS1(sim_reg, Cn2, D_0, d_0, dz, lambda)
%phase_screen_MVKS generate a phase screen built from the Modified Von Karman Spectrum (MVKS)
%   NOTE
%   Assuming uniform sampling and presents reflections on the boundaries
%
%   PARAMETERS
%   - sim_reg : SimulationRegion
%       SimulationRegion object containing the information about the simulation window
%   - Cn2 : float
%   - D_0 : float
%       large scale length
%   - d_0 : float
%       small scale length
%   - dz : float
%       thickness of the phase screen
%   - lambda : float
%       wavelength
%   
%   RETURNS
%   - phase_screen : complex
%       random phase screen generated from the MVKS

    arguments
        sim_reg SimulationRegion;
        Cn2 {mustBePositive};
        D_0 {mustBePositive};
        d_0 {mustBePositive};
        dz {mustBeNumeric};
        lambda {mustBePositive} = 1.55e-6;
    end

    % Define the sampling grid over the spatial frequency domain
    % Sampling inteval size
    dx = sim_reg.L_x/sim_reg.N_x;
    dy = sim_reg.L_x/sim_reg.N_y;
    fx = -1/(2*dx):1/sim_reg.L_x:(1/(2*dx)-1/sim_reg.L_x);
    fy = -1/(2*dy):1/sim_reg.L_y:(1/(2*dy)-1/sim_reg.L_y); 
    [Fx,Fy] = meshgrid(fx,fy);
    Kx = 2*pi*Fx;
    Ky = 2*pi*Fy;

    % Compute the MVKS
    k_m=5.92/d_0; % equivalent small scale wavenumber
    k_0=2*pi/D_0; % equivalent large scale wavenumber
    phi_n=0.033*Cn2./(Kx.^2+Ky.^2+ k_0^2).^(11/6).*exp(-(Kx.^2+Ky.^2).^2/k_m^2); % refractive index PSD according to the MVK model (eq. 11 [1])
    phi=2*pi*(2*pi/lambda)^2*dz*phi_n; % phase PSD (eq. 24 [1])

    % Generate a random seed
    seed=zeros(sim_reg.N_x, sim_reg.N_y)+3.14159265359i*ones(sim_reg.N_x, sim_reg.N_y);
    sigma_phi=(2*pi/(sim_reg.L_x))^2*phi;

    % Compute the phase screen by filtering the random seed through the MVKS spectrum 
    phase_screen=sim_reg.N_x^2*real(fftshift((ifft2(ifftshift(seed.*sqrt(sigma_phi))))));
end