function [r0] = r0(Cn2, L, k)
%r0 Fried parameter
%   Compute the fried parameter in the weak fluctuactions theory [1]
%   
%   PARAMETERS
%   - Cn2 : float
%   - L : float
%       channel length
%   - k : float
%       wavenumber
%   
%   RETURNS
%   - r0 : float
%       Fried parameter
%
%   REFERENCES
%   [1] L. C. Andrews and R. L. Phillips, Laser Beam Propagation through Random Media, 2nd ed. SPIE, 2005. doi: 10.1117/3.626196

    r0 = (0.423*k^2*Cn2*L).^(-3/5);
end

