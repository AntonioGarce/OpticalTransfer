    function [DOC] = DOC(U, r1, r2)
%DOC Degree of coherence
%   Compute the degree of coherence between the points two points, r1 and r2 starting from the mutual coherence
%   function as shown in [1]
%
%   PARAMETERS
%   - U : complex
%       3rd rank tensor containing different realizations of the complex electric
%       field. The last index identifies each realization.
%   - r1 : float
%       coordinates of the first point.
%   - r2 : float
%       coordinates of the second point.
%
%   RETURNS
%   - DOC : float
%       degree of coherence
%
%   REFERENCES
%   [1] L. C. Andrews and R. L. Phillips, Laser Beam Propagation through Random Media, 2nd ed. SPIE, 2005. doi: 10.1117/3.626196.

    arguments
        U (:, :, :);
        r1 (1, 2);
        r2 (1, 2);
    end

    DOC = abs(MCF(U, r1, r2))/sqrt(MCF(U, r1, r1)*MCF(U, r2, r2));
end