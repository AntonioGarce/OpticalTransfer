classdef SimulationRegion
    %SimulationRegion Provide information about the simulation region
    %   This object provide all the information about the simulation
    %   region, such as the window dimension, grid resolution, number of
    %   samples, etc.

    properties
        L_x;    % x-width of the simulation windows
        L_y;    % y-width of the simulation windows
        N_x;    % sample number along the horizontal axis
        N_y;    % sample number along the vertical axis
    end
    properties (Dependent, SetAccess = private)
        x;  % simulation window horizontal axis
        y;  % simulation window vertical axis
        X;  % meshgrid for the horizontal axis
        Y;  % meshgrid for the vertical axis
    end

    methods
        function obj = SimulationRegion(L_x, L_y, N_x, N_y)
            %SimulationRegion initialize the SimulationRegion object
            %   PARAMETERS
            %   - L_x : float
            %       x-width of the simulation windows
            %   - L_y : float
            %       y-width of the simulation windows
            %   - N_x : float
            %       sample number along the horizontal axis
            %   - N_y : float
            %       sample number along the vertical axis
            %
            %   RETURNS
            %       obj : SimulationRegion
            %           new instance of SimulationRegion

            obj.L_x = L_x;
            obj.L_y = L_y;
            obj.N_x = N_x;
            obj.N_y = N_y;
        end

        function value = get.x(obj)
        % getter method for the x property
            value = linspace(-obj.L_x/2, obj.L_x/2, obj.N_x);
        end

        function value = get.y(obj)
        % getter method for the y property
            value = linspace(-obj.L_y/2, obj.L_y/2, obj.N_y);
        end

        function value = get.X(obj)
        % getter method for the X property
            [value, ~] = meshgrid(obj.x, obj.y);
        end

        function value = get.Y(obj)
        % getter method for the Y property
            [~, value] = meshgrid(obj.x, obj.y);
        end
    end
end