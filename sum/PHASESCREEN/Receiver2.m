classdef Receiver2 < Receiver
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        N = 10; % Number of apertures
        center = [1, 1]; % center of the receiver
        w0 {mustBePositive} = 3e-3; % Aperture spot size
        is_symmetric logical = true;
        disp; % coordinates of the displacements
    end
    properties (Dependent)
        U;
    end

    methods
        function obj = Receiver2(simulation_region)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj@Receiver(simulation_region);
%             obj.disp = disp;
        end
    end
    methods
        function U = get.U(obj)
            if obj.is_symmetric
                U = exp(-((obj.sim_reg.X - obj.disp(1, 1)).^2 + (obj.sim_reg.Y - obj.disp(1, 2)).^2)/(obj.w0^2));
                U = U + exp(-((obj.sim_reg.X + obj.disp(1, 1)).^2 + (obj.sim_reg.Y + obj.disp(1, 2)).^2)/(obj.w0^2));
            else    
                U = exp(-((obj.sim_reg.X - obj.disp(1, 1)).^2 + (obj.sim_reg.Y - obj.disp(1, 2)).^2)/(obj.w0^2));
                U = U + exp(-((obj.sim_reg.X - obj.disp(2, 1)).^2 + (obj.sim_reg.Y - obj.disp(2, 2)).^2)/(obj.w0^2));
            end
            % normalize the receiver field
            U = U./sqrt(int2(abs(U).^2, obj.sim_reg));
        end
    end
end