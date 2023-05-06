function u_out= free_prop(u_in, sim_reg, z, lambda)
%Fresnel propagation using the Transfer function
% Based on Computational Fourier Optics by Voelz 
%
% Assuming uniform sampling and presents reflections on the boundaries
% 
%   PARAMETERS
%   - u_in : complex
%       complex amplitude of the beam at the object plane
%   - sim_reg : SimulationRegion
%       SimulationRegion object containing the informations on the simulation window.
%   - L : float
%       propagation length
%   - lambda : float
%       wavelength
%   RETURNS
%   - u_out : complex
%       complex aplitude of the beam at the image plane
%       
%   REFERENCES
%   [1] D. G. Voelz, Computational Fourier Optics: A MATLAB Tutorial. SPIE, 2011. doi: 10.1117/3.858456.

arguments
    u_in (:, :) {mustBeNumeric};
    sim_reg SimulationRegion;
    z {mustBeNumeric};
    lambda {mustBePositive} = 1.55e-6;
end

% Define the sampling grid over the spatial frequency domain
% Sampling inteval size
dx=sim_reg.L_x/sim_reg.N_x;
dy=sim_reg.L_x/sim_reg.N_y;
fx=-1/(2*dx):1/sim_reg.L_x:(1/(2*dx)-1/sim_reg.L_x);
fy=-1/(2*dy):1/sim_reg.L_y:(1/(2*dy)-1/sim_reg.L_y); 
[FX,FY]=meshgrid(fx,fy);

H=exp(-1j*pi*lambda*z*(FX.^2+FY.^2)); % Free propagation transfer function
H=fftshift(H);

U_in =fft2(fftshift(u_in)); % Angular spectrum of the input field
U_out=H.*U_in;
u_out=ifftshift(ifft2(U_out));  

end
