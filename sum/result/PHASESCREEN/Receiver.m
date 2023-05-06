classdef (Abstract) Receiver
    %Receiver Abstract class implementing a receiver
    %   

    properties (Abstract)
        N {mustBeInteger}; % Number of apertures
        center (1, 2) {mustBeReal}; % center of the receiver
    end
    properties (SetAccess = private)
        sim_reg SimulationRegion;
    end
    properties (Abstract, Dependent)
        U
    end

    methods
        function obj = Receiver(simulation_region)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.sim_reg = simulation_region;
        end
    end
end