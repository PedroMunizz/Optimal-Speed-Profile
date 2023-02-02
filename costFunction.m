function [J, P] = costFunction(x,d,theta)
    
    % J - Fuel consumption [l]    
    % x - Speed [m/s]
    % d - Distance of sector [m]
    % theta - Slope [rad]
    % vRol - Nominal Speed [m/s]
    
    % Parameters (Wang 2017)
       % Rolling (Rakha 2001)
       [nn,~] = size(x);
       Cr     = 1.25;                                     % Track roughness parameter
       c1     = 0.0328;                                   % Tire parameter 1
       c2     = 4.575;                                    % Tire parameter 2
       CFr    = ((ones(nn,1)*c2)+(3.6*c1*x))*(Cr/1000);   % Rolling resistance coefficient
        
	   % Aero
	   rho = 1.2256;                                       % Air density [kg/m3]
	   Cd  = 0.7;                                          % Drag coefficient
	   A   = 8.9;                                          % Frontal area [m2]
	   C   = 0.5*rho*Cd*A;                                 % Lumped air drag coefficient

	   % Inertial
	   m = 16000;                                          % Vehicle mass [kg]
        
       % Gravity
	   g = 9.81;                                           % [m/s2]
       efPower = 0.915;                                    % Powertrain efficiency
       
	   % Fuel consumption model (Wang 2017)
	   alpha0 = 1.87E-01;
	   alpha1 = 5.23E-02;
	   alpha2 = 3.47E-04;

    % Resistive force
    R = C*x.^2 + m*g*sin(theta) + m*g*CFr;

    % Power
    P = (R.*x)/(1000*efPower);
    
    % Fuel rate
    	% Preallocating
    	FC = zeros(length(x),1);

        for i=1:length(x)
            if P(i)>0
                FC(i) = alpha0 + alpha1*P(i) + alpha2*P(i)^2;
            else
                FC(i) = alpha0;
            end
        end

    % Cost
	J = (sum(FC.*d./x))/1000;         %kg/s
    
end
