function c = overlap(U_2, U_1, sim_reg)
% overlap mode overlap between two optical fields under scalar approximation
%   
%   PARAMETERS
%   - U_2 : complex
%       complex amplitude of the optical field to project onto
%   - U_1 : complex
%       complex amplitude of the optical field projected onto U_1
%   - sim_reg : SimulationRegion
%       SimulationRegion object containing the information about the simulation window
%   
%   RETURNS
%   - c : complex
%       Coupling coefficient

    arguments
        U_2;
        U_1;
        sim_reg SimulationRegion;
    end
    
    % WARNING: this may work only for normalized fields
    c = int2(conj(U_2).*U_1, sim_reg);
end