function [c, ceq] = nonLinConst(x,v,d,theta,vNom)
 
    % x - Speed [m/s]
    % v - Nominal Speed [m/s]
    % d - Distance of equally spaced sectors [m]
    
    % Sector speeds
    distT = sum(d);         % Total length of the stretch

    % Time and acceleration matrix for ramps
    [nn,~] = size(x);   
       
       % Rolling (Rakha 2001)
       Cr     = 1.25;                                      % Track roughness parameter
       c1     = 0.0328;                                    % Tire parameter 1
       c2     = 4.575;                                     % Tire parameter 2
       CFr    = ((ones(nn,1)*c2)+(3.6*c1*x))*(Cr/1000);    % Rolling resistance coefficient
        
	   % Aero
	   rho = 1.2256;                                       % Air density [kg/m3]
	   Cd  = 0.7;                                          % Drag coefficient
	   A   = 8.9;                                          % Frontal area [m2]
	   C   = 0.5*rho*Cd*A;                                 % Lumped air drag coefficient

	   % Inertial
	   m      = 16000;                                     % Vehicle mass [kg]
        
       % Gravity
	   g = 9.81;                                           % [m/s2]
       efPower = 0.915;                                    % Powertrain efficiency [-]

    % Resistive force
     R = C*x.^2 + m*g*sin(theta) + m*g*CFr;

    % Power
     P1 = (R.*x) / (1000*efPower);  %Power
     P1_length = length(P1);
     Pmax = ones(P1_length,1)*150;
   
    % Inequality constraints (c <= 0)
    c = P1-Pmax;

    % Equality constraints (ceq = 0)
    ceq = [sum(d./x) - distT/v; x(1) - vNom; x(end) - vNom]; % Total time and initial v identical
    
end



