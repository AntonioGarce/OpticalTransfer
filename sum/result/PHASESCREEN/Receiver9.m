classdef Receiver9 < Receiver
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        N = 2; % Number of apertures
        center = [0, 0]; % center of the receiver
        w0 {mustBePositive} = 3e-3; % Aperture spot size
        is_symmetric logical = true;
        disp; % coordinates of the displacements
    end
    properties (Dependent)
        U;
    end

    methods
        function obj = Receiver9(simulation_region, disp)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj@Receiver(simulation_region);
            obj.disp = disp;
        end
    end
    methods
        function U = get.U(obj)
            if obj.is_symmetric
                 U = exp(-((obj.sim_reg.X - obj.disp(1, 1)).^2 + (obj.sim_reg.Y - obj.disp(1, 2)).^2)/(obj.w0^2));
%                     U = U + exp(-((obj.sim_reg.X + obj.disp(1, 1)).^2 + (obj.sim_reg.Y + obj.disp(1, 2)).^2)/(obj.w0^2));
                disp_v = obj.disp(1, 1)+1j*obj.disp(1, 2);   
                for i=1:8
                    disp_v = disp_v * (cos(2*pi/9)+1j*sin(2*pi/9));
                    xc = real(disp_v);
                    yc = imag(disp_v);
                     U = U + exp(-((obj.sim_reg.X - xc).^2 + (obj.sim_reg.Y - yc).^2)/(obj.w0^2));
                end
            else    
                for k = 1:9
                     U = exp(-((obj.sim_reg.X - obj.disp(k, 1)).^2 + (obj.sim_reg.Y - obj.disp(k, 2)).^2)/(obj.w0^2));
                end
%                 U = exp(-((obj.sim_reg.X - obj.disp(1, 1)).^2 + (obj.sim_reg.Y - obj.disp(1, 2)).^2)/(obj.w0^2));
%                 U = U + exp(-((obj.sim_reg.X - obj.disp(2, 1)).^2 + (obj.sim_reg.Y - obj.disp(2, 2)).^2)/(obj.w0^2));
            end
            % normalize the receiver field
            U = U./sqrt(int2(abs(U).^2, obj.sim_reg));
        end
    end
end