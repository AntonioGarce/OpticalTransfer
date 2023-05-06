classdef (Abstract) Receiver
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        N mustBeInteger = 1; % Number of apertures
        center (1, 2) mustBeReal = [0, 0]; % center of the receiver
        s_x = 1;
        s_y = 1;
        w0 mustBePositive = 0; % Aperture spot size
        piano;
        is_circular boolean;
        is_equispaced boolean;
    end
    properties (SetAccess = private)
        sim_reg SimulationRegion;
    end
    properties (Abstract, Dependent)
        U
    end

    methods (abstract)
        function obj = Receiver(simulation_region)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.sim_reg = simulation_region;
        end

        U = get.U(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            U = zeros(size(obj.sim_reg.X));
                switch N
                    case 1
                        U = U + exp((-(obj.sim_reg.X-obj.center(1,1)-sx).^2-(obj.sim_reg.Y-obj.center(1,2)).^2)/(obj.w0^2));
                        
                    case 2
                        
                        if symmetric==1
                            k=1;
                            for i=1:N
                                if i==2
                                    k=-1;
                                end
                                U = U + exp((-(X-k*sx/2).^2-(Y-k*sy/2).^2)/(w0^2));
                            end
                            U = U./sqrt(trapz(assey,trapz(assex,(abs(U)).^2,2)));
                        elseif symmetric==0
                            for i=1:N
                                U = U + exp((-(X-center(1,1)-(i-1)*sx).^2-(Y-center(1,2)).^2)/(w0^2));
                            end
                            U = U./sqrt(trapz(assey,trapz(assex,(abs(U)).^2,2)));
                        end
                        
                        
                    case 4
                        k = [1,1; -1,-1; 1,-1; -1,1];
                        for i=1:N
                            U = U + exp((-(X-k(i,1)*sx*sqrt(2)/2).^2-(Y-k(i,2)*sx*sqrt(2)/2).^2)/(w0^2));
                        end
                        U = U./sqrt(trapz(assey,trapz(assex,(abs(U)).^2,2)));
                        
                        
                    case 5
                        
                        k = [0,0; 1,1; -1,-1; 1,-1; -1,1];
                        for i=1:N
                            U = U + exp((-(X-k(i,1)*sx*sqrt(2)/2).^2-(Y-k(i,2)*sx*sqrt(2)/2).^2)/(w0^2));
                        end
                        U = U./sqrt(trapz(assey,trapz(assex,(abs(U)).^2,2)));
                        
                        
                    case 8   % di default i pallini sono distanziati di 4e-2 (direzione radiale)
                        p=8e-2;
                        k = [ 0,0; -0.04*p*1.25,0.39*p*1.25; 0.28*p*1.25,0.28*p*1.25; 0.39*p*1.25,-0.04*p*1.25; 0.21*p*1.25,-0.33*p*1.25; -0.132112*p*1.25,-0.377553*p*1.25; -0.37*p*1.25,-0.13*p*1.25; -0.33*p*1.25,0.21*p*1.25;];   %7
                        for i=1:N
                            U = U + exp((-(X-k(i,1)*sx).^2-(Y-k(i,2)*sx).^2)/(w0^2));
                        end
                        U = U./sqrt(trapz(assey,trapz(assex,(abs(U)).^2,2)));
                
                        
                    case 9
                        if cerchio == 1
                            k = [0,0; 1,0; 0,1; -1,0; 0,-1; 1/sqrt(2),1/sqrt(2); -1/sqrt(2),-1/sqrt(2); 1/sqrt(2),-1/sqrt(2); -1/sqrt(2),1/sqrt(2)];
                            for i=1:N
                                U = U + exp((-(X-k(i,1)*sx).^2-(Y-k(i,2)*sx).^2)/(w0^2));
                            end
                            U = U./sqrt(trapz(assey,trapz(assex,(abs(U)).^2,2)));
                        elseif cerchio == 0
                            k = [0,0; 0,2; 0,-2; -2,0; 2,0; 1/sqrt(2),1/sqrt(2); -1/sqrt(2),-1/sqrt(2); 1/sqrt(2),-1/sqrt(2); -1/sqrt(2),1/sqrt(2)];
                            for i=1:N
                                U = U + exp((-(X-k(i,1)*sx).^2-(Y-k(i,2)*sx).^2)/(w0^2));
                            end
                            U = U./sqrt(trapz(assey,trapz(assex,(abs(U)).^2,2)));
                        end
                        
                        
                    case 16  % di default i pallini sono distanziati di 4e-2 (direzione radiale) tra loro
                        if equispaziato == 1
                            p=8e-2;
                            k = [ 0,0; p,0; 0,p; -p,0; 0,-p; 1/sqrt(2)*p,1/sqrt(2)*p; -1/sqrt(2)*p,-1/sqrt(2)*p; 1/sqrt(2)*p,-1/sqrt(2)*p; -1/sqrt(2)*p,1/sqrt(2)*p;  %9
                                -0.04*p*1.25,0.39*p*1.25; 0.28*p*1.25,0.28*p*1.25; 0.39*p*1.25,-0.04*p*1.25; 0.21*p*1.25,-0.33*p*1.25; -0.132112*p*1.25,-0.377553*p*1.25; -0.37*p*1.25,-0.13*p*1.25; -0.33*p*1.25,0.21*p*1.25;];   %7
                            for i=1:N
                                U = U + exp((-(X-k(i,1)*sx).^2-(Y-k(i,2)*sx).^2)/(w0^2));
                            end
                            
                            U = U./sqrt(trapz(assey,trapz(assex,(abs(U)).^2,2)));
                            
                        elseif equispaziato == 0
                            p=8e-2;
                            k = [ 0,0; p,0; 0,p; -p,0; 0,-p; 1/sqrt(2)*p,1/sqrt(2)*p; -1/sqrt(2)*p,-1/sqrt(2)*p; 1/sqrt(2)*p,-1/sqrt(2)*p; -1/sqrt(2)*p,1/sqrt(2)*p;  %9
                                -0.04*p*2/3*1.25,0.39*p*2/3*1.25; 0.28*p*2/3*1.25,0.28*p*2/3*1.25; 0.39*p*2/3*1.25,-0.04*p*2/3*1.25; 0.21*p*2/3*1.25,-0.33*p*2/3*1.25; -0.132112*p*2/3*1.25,-0.377553*p*2/3*1.25; -0.37*p*2/3*1.25,-0.13*p*2/3*1.25; -0.33*p*2/3*1.25,0.21*p*2/3*1.25;];   %7
                            for i=1:N
                                U = U + exp((-(X-k(i,1)*sx).^2-(Y-k(i,2)*sx).^2)/(w0^2));
                            end
                            
                            U = U./sqrt(trapz(assey,trapz(assex,(abs(U)).^2,2)));
                            
                        end
                        
                        
                    case 17
                        k = [0,0; 1,0; 0,1; -1,0; 0,-1; 1/sqrt(2),1/sqrt(2); -1/sqrt(2),-1/sqrt(2); 1/sqrt(2),-1/sqrt(2); -1/sqrt(2),1/sqrt(2); -2,0; 2,0; 0,2; 0,-2;  2/sqrt(2),2/sqrt(2); -2/sqrt(2),-2/sqrt(2); 2/sqrt(2),-2/sqrt(2); -2/sqrt(2),2/sqrt(2)];
                        
                        for i=1:N
                            U = U + exp((-(X-k(i,1)*sx).^2-(Y-k(i,2)*sx).^2)/(w0^2));
                        end
                end
                % normalize the receiver field
                U = U./sqrt(int2(U, obj.sim_reg));
        end
    end
end