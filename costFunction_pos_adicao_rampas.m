function [J, time_add, d_f, J_constant, J_ramps, P1, P2, d2] = costFunction_pos_adicao_rampas(x,accel,d,delta,theta)
    
    % J     - Fuel consumption [l]    
    % x     - Optimal speed profile [m/s]
    % d     - Distance of sector [m]
    % acel  - Constant acceleration for ramps [m/s^2]
    % theta - Slope [rad]
    % delta - Initial step for the integration
    % vRol - Nominal Speed [m/s]
    
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
       lambda = 0.1;                                       % Rotational masses mass factor [-]
        
       % Gravity
	   g = 9.81;                                           % [m/s2]
       efPower = 0.915;                                    % Powertrain efficiency [-]
       
	   % Fuel consumption model (Wang 2017)
	   alpha0 = 1.87E-01;
	   alpha1 = 5.23E-02;
	   alpha2 = 3.47E-04;

    % Resistive force
     R = C*x.^2 + m*g*sin(theta) + m*g*CFr;

    % Power
     P1 = (R.*x) / (1000*efPower);  %Power using the first vel. of the interval
    
    % Fuel rate (constant speed)
    	% Preallocating
    	FC = zeros(length(x),1);

        for i=1:length(x)
            if P1(i)>0
                FC(i) = alpha0 + alpha1*P1(i) + alpha2*P1(i)^2;
            else
                FC(i) = alpha0;
            end
        end
        
    % Fixing d matrix for correct consumption of J_constant
    d_ramps = zeros((nn),1);                               % Preallocating distances for ramps consumption
    
    for i=2:length(d_ramps)                                % Time matrix for ramps consumption
        d_ramps(i)= (x(i)^2 - x(i-1)^2)/(2*a_mtx(i-1));    % Distance = (v2^2 - v1^2) / (2*a)
    end
    
    d_f = d - d_ramps;
        
     % Fuel rate (acceleration ramps)
    	% Preallocating
    	Ramps_J = zeros((nn-1),1);
        hold = 0;

        for i=1:length(t_ramps)
            l = ceil( abs(x(i+1)-x(i)) / delta );
            delta_new = (x(i+1)-x(i)) / l;
            iterat_t = delta_new / a_mtx(i);
            v_atual = x(i);
            hold = hold + d_f(i);                          %Hold the distance value until the other loop
            for j=1:l
                R = C*v_atual^2 + m*g*sin(theta(i+1)) + m*g*(Cr/1000)*(c2+(3.6*c1*v_atual));
                Piterat = ((R + m*(1+lambda)*a_mtx(i))*v_atual) / (1000*efPower);
                if abs(x(i+1)-x(i)) < 0.01
                    P2(l*(i-1)+j) = (R*v_atual) / (1000*efPower);
                else
                    P2(l*(i-1)+j) = Piterat;
                end
                d2(l*(i-1)+j) = hold + v_atual*iterat_t;
                hold = d2(l*(i-1)+j);

                if Piterat>0
                    iterat_FC = alpha0 + alpha1*Piterat + alpha2*Piterat^2;
                    iterat_J = iterat_FC * iterat_t;
                    Ramps_J(i) = Ramps_J(i) + iterat_J;
                else
                    iterat_FC = alpha0;
                    iterat_J = iterat_FC * iterat_t;
                    Ramps_J(i) = Ramps_J(i) + iterat_J;
                end
                v_atual = v_atual + delta_new;
            end
        end
    
    % Cost
 	J_constant = (sum(FC.*d_f./x)) / 1000;                 %kg/s
	J_ramps = (sum(Ramps_J)) / 1000;                       %kg/s
    time_add = t_ramps;
    J = J_constant + J_ramps;   
end