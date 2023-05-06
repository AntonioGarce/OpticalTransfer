function [Cn2] = Cn(sigma_ritov, k, L)
% Cn compute the C_n^2> parameter
%   Compute the C_n^2 parameter from the Rytov variance as in [1]
%
%   PARAMETERS
%   - sigma_ritov : float
%       Rytov variance
%   - k : float
%       wavenumber of the laser beam
%   - L : float
%       length of the link
%
%   RETURNS
%   - Cn2 : float
%       C_n^2 parameter
%   
%   REFERENCES
%   [1] L. C. Andrews and R. L. Phillips, Laser Beam Propagation through Random Media, 2nd ed. SPIE, 2005. doi: 10.1117/3.626196.

    arguments
        sigma_ritov {mustBePositive};
        k {mustBePositive};
        L {mustBePositive};
    end

    Cn2 = sigma_ritov./(1.23*k^(7/6)*L^(11/6));
end

