function [sigma] = scint_idx(img)
% scint_idx scintillation index
%   Compute the scintillation indx from a collection of optical field realizations
%   
%   PARAMETERS
%   - U : complex
%       3rd rank tensor containing different realizations of the complex electric
%       field. The last index identifies each realization.
%
%   RETURNS
%   -sigma : complex
%       scintillation index computed point by point

    sigma = mean(img.^2,3)./((mean(img,3)).^2) - 1; 
end