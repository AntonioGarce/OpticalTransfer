function [U_out, phase_screen_used] = turbulent_prop_new(U_in, sim_reg, L, Cn2, D0, d0, n_screen, phase_screen_old, corr_coeff,lamdba)
% turbulent_prop laser beam propagation in a turbulent medium 
%   Compute the propagation of a laser beam through a turbulent medium via the split-step method
%
%   PARAMETERS
%   - U_in : complex
%       complex amplitude of the electromagnetic field at the object plane
%   - sim_reg : SimulationRegion
%       SimulationRegion object containing the information about the simulation window
%   - L : float
%       propagation length
%   - Cn2 : float
%   - D_0 : float
%       large scale length
%   - d_0 : float
%       small scale length
%   - dz : float
%       thickness of the phase screen
%   - n_screen : uint
%       number of phase screen employed
%   - phase_screen_old : complex
%       phase screens from the previous iteration, necessary to generate temporal correlated turbulence
%   - corr_coeff : float
%       correlation coefficient
%   - lambda : float
%       wavelength
%
%   RETURNS
%   - U_out : complex
%       complex amplitude of the electromagnetic field at the object plane
%   - phase_screen_used : complex
%       phase screen used during the execution of the split step method

    arguments
        U_in (:, :) {mustBeNumeric}; % Input beam
        sim_reg SimulationRegion; % Simulation region
        L {mustBePositive}; % Length of the channel
        Cn2 {mustBePositive};
        D0 {mustBePositive}; % Outer scale of the MVK model
        d0 {mustBePositive}; % inner scale of the MVK model
        n_screen {mustBeInteger, mustBePositive}; % Number of phase screen used
        phase_screen_old;
        corr_coeff {mustBeLessThan(corr_coeff, 1), mustBePositive} = 0; % Correlation coefficient of the output phase screen with the input one
        lamdba {mustBePositive} = 1.55e-6;
    end
    
    % Generate phase screens
    dL = L/(n_screen + 1);
    phase_screen_used = zeros(sim_reg.N_x, sim_reg.N_y, n_screen);
    for i = 1:n_screen
        phase_screen_used(1:(sim_reg.N_y)/2, :, i) = zeros;
        phase_screen_used((sim_reg.N_y)/2+1:sim_reg.N_y, :, i) = pi;
    end
    % Correlate the phase screen with those used in a previous iteration
    phase_screen_used = sqrt(1-corr_coeff^2)*phase_screen_used + corr_coeff*phase_screen_old;
    
    U_out = U_in;
    for i=1:n_screen
        U_out = free_prop(U_out, sim_reg, dL); %propagation in a dz section
        U_out = U_out.*exp(1i*phase_screen_used(:, :, i));
    end
    % Propagation after the last phase screen
    U_out = free_prop(U_out, sim_reg, dL);
end