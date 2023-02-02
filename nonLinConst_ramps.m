function [c, ceq] = nonLinConst_ramps(x,v,d,accel,theta,vNom)
 
    % x - Speed [m/s]
    % v - Nominal Speed [m/s]
    % d - Distance of equally spaced sectors [m]
    
    % Sector speeds
    distT = sum(d);                                       % Total length of the stretch
    
    % Time and acceleration matrix for ramps
    [nn,~] = size(x);
    t_ramps = zeros((nn-1),1);                            % Preallocating times for ramps consumption
    a_mtx = zeros(nn,1);                                  % Preallocating acceleration matrix
       
       for i=1:length(t_ramps)                            % Time matrix for ramps consumption
           if x(i+1)-x(i) <0
               t_ramps(i) = (x(i+1)-x(i))/(-accel);
               a_mtx(i) = -1;
           else
               t_ramps(i) = (x(i+1)-x(i))/accel;
               a_mtx(i) = 1;
           end
       end
    
    a_mtx = a_mtx*accel;                                  % Acceleration matrix (pos & neg numbers)
       
    % Distance matrix for ramps
    d_ramps = zeros((nn),1);                              % Preallocating distances for ramps consumption
        for i=2:length(d_ramps)                           % Time matrix for ramps consumption
           d_ramps(i)= (x(i)^2 - x(i-1)^2)/(2*a_mtx(i-1));% Distance = (v2^2 - v1^2) / (2*a)
        end
    
    %Total time spent on ramps
    d_f = d - d_ramps;
    t_vcte = sum(d_f./x);
      
    % Parameters (Wang 2017)
       [nn,~] = size(x);
       t_ramps = zeros((nn-1),1);                          % Preallocating times for ramps consumption
       a_mtx = zeros(nn,1);                                % Preallocating acceleration matrix
       
       for i=1:length(t_ramps)                             % Time matrix for ramps consumption
           if x(i+1)-x(i) <0
               t_ramps(i) = (x(i+1)-x(i))/(-accel);
               a_mtx(i) = -1;              
           else
               t_ramps(i) = (x(i+1)-x(i))/accel;
               a_mtx(i) = 1;
           end
       end
      
       a_mtx = a_mtx*accel;                                 % Acceleration matrix (pos & neg numbers)
       
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
     d_f_length = length(d_f);
     dmin = ones(d_f_length,1)*100;
    
    % Inequality constraints (c <= 0)
    c = [(t_vcte + sum(t_ramps))-(distT/v); P1-Pmax; dmin-d_f];
    
    % Equality constraints (ceq = 0)
    ceq = [x(1) - vNom; x(end) - vNom];
    
end


