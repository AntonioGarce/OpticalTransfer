function [MCF] = MCF(U, r1, r2)
%MCF Mutual coherence function
%   Compute the mutual coherence function o a stochastic optical field U between two points r1 and r2
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
%   - MCF : float
%       mutual coherence function

    arguments
        U (:, :, :);
        r1 (1, 2);
        r2 (1, 2);
    end
    
        
    MCF = mean(U(r1(1), r1(2), :).*conj(U(r2(1), r2(2), :)));
end
