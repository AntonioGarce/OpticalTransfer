classdef Receiver1 < Receiver
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        N = 1; % Number of apertures
        center = [0, 0]; % center of the receiver
        w0 {mustBePositive} = 3e-8; % Aperture spot size
    end
    properties (Dependent)
        U;
    end

    methods
        function obj = Receiver1(simulation_region)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj@Receiver(simulation_region);
        end
    end
    methods
        function U = get.U(obj)
            U = exp(-((obj.sim_reg.X - obj.center(1,1)).^2 + (obj.sim_reg.Y - obj.center(1,2)).^2)/(obj.w0^2));

            % normalize the receiver field
            U = U./sqrt(int2(U, obj.sim_reg));
        end
    end
end