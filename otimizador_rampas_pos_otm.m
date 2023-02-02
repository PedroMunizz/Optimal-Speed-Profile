%% Problema de controle ótimo

%% Defining optimization scenario

dTotal  = 11000;                     % Total distance of the stretch - Trajeto A (Lombada) [m]
%dTotal  = 8000;                      % Total distance of the stretch - Trajeto B (Sigmoide) [m]
%dTotal  = 16000;                     % Total distance of the stretch - Trajeto C Washington Luis [m]
%dTotal  = 4000;                      % Total distance of the stretch - Trajeto D Presidente Dutra [m]
%dTotal  = 26000;                     % Total distance of the stretch - Trajeto E Castello Branco [m]
vNom    = 55/3.6;                    % Nominal speed [m/s]
dSector = 500;                       % Distance of each sector [m]
Nsec    = round(dTotal/dSector);     % Number of sectors
accel   = 0.15;                      % Constant acceleration for ramps [m/s^2]
delta   = 0.005;                     % Initial delta (step) for the integration between two velocities
c1 = 60; c2 = 6000; c3 = 1150;        % Coefficients for Trajeto A - Lombada
csig1 = 50; csig2 = 300;             % Coefficients for Trajeto B - Sigmoide
Ttotal = dTotal/vNom;                % Total period [s]

% Generating equally spaced sectors
sect = ones(Nsec,1)*dTotal/Nsec;
    
% Defining the average constant slope of each sector
    % Trajeto A - Lombada (usar dTotal = 11000 m)
    acum = 0;
    altit = zeros(Nsec+1,1);
    theta = zeros(Nsec,1);
    for i=1:Nsec
        acum = acum + dSector;
        altit(i+1) = c1*exp(-((acum-c2)^2/(2*c3^2)));
    end

    for i=1:Nsec
        theta(i) = asin((altit (i+1) - altit(i)) / dSector);
    end

    % Trajeto B - Sigmoide (usar dTotal = 8000 m)
%     acum = -dTotal/2;
%     altit = zeros(Nsec+1,1);
%     theta = zeros(Nsec,1);
%     for i=1:Nsec
%         acum = acum + dSector;
%         altit(i+1) = csig1 / (1+exp(-acum/csig2));
%     end
%     
%     for i=1:Nsec
%         theta(i) = asin((altit (i+1) - altit(i)) / dSector);
%     end    

    % Trajeto C - Washington Luis (usar dTotal = 16000 m)
    %theta = [-0.014;-0.04;-0.056;-0.004;-0.01;-0.012;-0.022;-0.002;-0.002;-0.008;0.018;0.038;0.012;0.012;0.016;-0.028;-0.026;-0.012;-0.018;-0.018;0.006;0.006;0.01;0.004;-0.006;-0.01;0.042;0;-0.002;0.024;-0.038;0.012];

    % Trajeto D - Presidente Dutra (usar dTotal = 4000 m)
    %theta = [0.0196; -0.0252; -0.0126; 0.0212; 0.0104; -0.0344; -0.0234; 0.0352];

    % Trajeto E - Castello Branco (usar dTotal = 26000 m)
%     theta = [-0.015;0.018;0.024;0.0168;-0.005;-0.021;-0.010;-0.006;-0.008;-0.001;-0.010;-0.031;-0.008;0.011;0.024;-0.011;-0.014;-0.008;-0.039;-0.011;-0.001;0.011;0.025;0.009;0.008;0.006;0.004;0.040;-0.007;-0.028;0.031;0.009;-0.007;0.010;-0.011;0.002;0.029;0;
%         -0.010;-0.011;-0.015;-0.016;0.009;-0.004;0.008;0.012;0.001;-0.006;-0.012;0.004;0.003;-0.007];

% Defining speed limits for the entire stretch
vMinVal = 45/3.6;                    % Minimum speed [m/s]
vMaxVal = 90/3.6;                    % Maximum speed [m/s]

vMin = ones(Nsec,1)*vMinVal;
vMax = ones(Nsec,1)*vMaxVal;

%% Running optimization using fmincon
x0 = ones(Nsec,1)*vNom;

% Array of lower and upper bounds
lb = vMin;
ub = vMax;

% fmincon options
options = optimoptions(@fmincon,'OptimalityTolerance',1.0000e-04,'MaxFunctionEvaluations',80000,'MaxIterations',8000);

% Running fmincon
[x_step,c_step] = fmincon(@(x) costFunction(x,sect,theta),x0,[],[],[],[],lb,ub,@(x) nonLinConst(x,vNom,sect,theta,vNom),options);

%fmincon for 2 initial assumptions: x0 e x_step
[x_accel,c_accel] = fmincon(@(x) costFunction_pos_adicao_rampas(x,accel,sect,delta,theta),x0,[],[],[],[],lb,ub,@(x) nonLinConst_ramps(x,vNom,sect,accel,theta,vNom),options);
[x_accel2,c_accel2] = fmincon(@(x) costFunction_pos_adicao_rampas(x,accel,sect,delta,theta),x_step,[],[],[],[],lb,ub,@(x) nonLinConst_ramps(x,vNom,sect,accel,theta,vNom),options);
if c_accel < c_accel2
    initial_assump = 'x0';
    clear x_accel2
    clear c_accel2
else
    initial_assump = 'x_step';
    c_accel = c_accel2;
    x_accel = x_accel2;
    clear x_accel2
    clear c_accel2
end

%% Verification

% Nominal cost
[Jnom, Pnom] = costFunction(ones(Nsec,1)*vNom,sect,theta);
Jnom

% Optimal cost post accel.|decel. ramps addition
[~,timeadd_m,d_f,J_constant,J_ramps, P1, P2, d_ramps] = costFunction_pos_adicao_rampas(x_accel,accel,sect,delta,theta);

% Savings
savings = 100*(Jnom-c_accel)/Jnom;
 
% Figures
[timeNomVec,distNomVec,veloNomVec,timeOptSTEP,veloOptSTEP,distOptSTEP,timeOptACCEL,veloOptACCEL,distOptACCEL,slopeOptVector,POptVector] = results(x_step,x_accel,dTotal,sect,vNom,accel,theta,P1);
velOPTaccel_km = transpose(veloOptACCEL .* 3.6);
velOPTstep_km = veloOptSTEP .* 3.6;
x_step_km = x_step .*3.6;

figure
    str0 = ['Total distance: ', num2str(dTotal),' m'];
    annotation('textbox',[.134 .98 .0 .0],'String',str0,'FitBoxToText','on');
    str1 = ['Nominal consumption: ', num2str(Jnom), ' kg'];
    annotation('textbox',[.134 .915 .0 .0],'String',str1,'FitBoxToText','on');
    str2 = ['Optimal consumption (step): ', num2str(c_step), ' kg'];
    annotation('textbox',[.134 .85 .0 .0],'String',str2,'FitBoxToText','on');
    str3 = ['Optimal consumption (ramps): ', num2str(c_accel), ' kg'];
    annotation('textbox',[.134 .785 .0 .0],'String',str3,'FitBoxToText','on');
figure
subplot(3,1,1)
    hold on ; box on ; grid on    
    title ('Curva de potência')
    plot(distOptACCEL/1000,POptVector,'-o','MarkerSize',4)
    plot(d_ramps/1000,P2,'o','MarkerSize',4)
    plot([0 dTotal/1000],[150 150])
    ylabel('Power [kW]')
subplot(3,1,2)
    hold on ; box on ; grid on
    title('Perfil de velocidade')
    plot(distOptSTEP/1000,veloOptSTEP*3.6)    
    plot(distOptACCEL/1000,veloOptACCEL*3.6)
    plot(distNomVec/1000,veloNomVec*3.6)
    ylabel('Velocidade [km/h]')
    legend('Optimal','Optimal_R_a_m_p_s','Nominal','Location','Best')
subplot(3,1,3)
    hold on ; box on ; grid on
    title('Perfil de relevo')
    plot(distOptSTEP/1000,tan(slopeOptVector)*100)
    ylabel('Inclinação [%]')
    xlabel('Distância [km]')