function [double_integral] = int2(f_samples, domain)
%int2 perform a double integral over a specified domain
%   Perform a double integral over a specified domain using a trapezoidal interpolation
%   
%   PARAMETERS
%   - f_samples : complex
%       values of the integrand functions at the interpolation nodes specified in domain
%   - domain : SimulationRegion
%       SimulationRegion object defining the integration domain and the interpolation nodes coordinates
%   
%   RETURNS
%   - double_integral : complex
%       double integral of the function over the specified domain

    arguments
        f_samples (:, :)
        domain SimulationRegion
    end

    double_integral = trapz(domain.x, trapz(domain.y, f_samples));
end

