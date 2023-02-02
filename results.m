function [timeNomVec,distNomVec,veloNomVec,timeOptSTEP,veloOptSTEP,distOptSTEP,timeOptACCEL,veloOptACCEL,distOptACCEL,slopeOptVector,POptVector] = results(xOpt_step,xOpt_accel,distTotal,d,velNom,accel,slope,P)
	% Generating the graphics for N sectors
    % xOpt_step - optimal step speed profile
    % xOpt_accel - optimal speed profile with acceleration ramps
    
	% Assuming equal distances
	N = length(xOpt_step);          % Number of sectors
	distSector = distTotal/N;       % Distance of each sector for STEP speed profile

    % Time and acceleration matrix for ramps
    t_ramps = zeros((N-1),1);                            % Preallocating times for ramps consumption
    a_mtx = zeros(N,1);                                  % Preallocating acceleration matrix
    for i=1:length(t_ramps)                           % Time matrix for ramps consumption
        if xOpt_accel(i+1)-xOpt_accel(i) <0
            t_ramps(i) = (xOpt_accel(i+1)-xOpt_accel(i))/(-accel);
            a_mtx(i) = -1;
        else
            t_ramps(i) = (xOpt_accel(i+1)-xOpt_accel(i))/accel;
            a_mtx(i) = 1;
        end
    end
    
    a_mtx = a_mtx*accel;   
    
    % Distance matrix for ramps
    d_ramps = zeros((N),1);                                                 % Preallocating distances for ramps consumption 
    for i=2:length(d_ramps)                                                 % Time matrix for ramps consumption
        d_ramps(i)= (xOpt_accel(i)^2 - xOpt_accel(i-1)^2)/(2*a_mtx(i-1));   % Distance = (v2^2 - v1^2) / (2*a)
    end
    
    d_f = d - d_ramps;
       
	% Nominal trajectory
	timeNomVec = [0 distTotal/velNom;]; % Time array nominal case
	distNomVec = [0 distTotal];			% Distance array nominal case
	veloNomVec = [velNom velNom];       % Velocity array nominal case

	% Optimal trajectory
		% Position curve
			% Preallocating
			distOptVecSTEP     = zeros(N+1,1);
			timeOptVecSTEP     = zeros(N+1,1);

			for i = 2:N+1
				vIterationSTEP = xOpt_step(i-1);
				timeOptVecSTEP(i) = timeOptVecSTEP(i-1) + distSector/vIterationSTEP;
				distOptVecSTEP(i) = distOptVecSTEP(i-1) + distSector;
            end
            
	%Velocity curve, slope, gear and RPM
		%Preallocating
			slopeOptVector   = zeros(2*N,1);
            POptVector       = zeros(2*N,1);
%           iUsoOptVector    = zeros(2*N,1);
%           wEngineOptVector = zeros(2*N,1);
			veloOptSTEP      = zeros(2*N,1);
			timeOptSTEP      = zeros(2*N,1);
			timeOptSTEP(end) = timeOptVecSTEP(end);
			distOptSTEP      = zeros(2*N,1);
			distOptSTEP(end) = distOptVecSTEP(end);
            distOptACCEL     = zeros(2*N,1);
            timeOptACCEL      = zeros(2*N,1);

			for i = 1:2:2*N
				vIterationSTEP = xOpt_step((i+1)/2);
				veloOptSTEP(i:i+1) = [vIterationSTEP vIterationSTEP];
                
                vIterationACCEL = xOpt_accel((i+1)/2);
				veloOptACCEL(i:i+1) = [vIterationACCEL vIterationACCEL];

				slopeIteration = slope((i+1)/2);
				slopeOptVector(i:i+1) = [slopeIteration slopeIteration];
                
                PIteration = P((i+1)/2);
				POptVector(i:i+1) = [PIteration PIteration];
%                 
%               iUsoIteration = ratio((i+1)/2);
% 				iUsoOptVector(i:i+1) = [iUsoIteration iUsoIteration];
%                 
%               wEngineIteration = wspeed((i+1)/2);
% 				wEngineOptVector(i:i+1) = [wEngineIteration wEngineIteration];             
            end
            
        %Creating the vector for the plots (duplicate the values)   
			for i = 2:2:2*N-1
                timeOptSTEP(i:i+1) = [timeOptVecSTEP(i/2+1) timeOptVecSTEP(i/2+1)];
				distOptSTEP(i:i+1) = [distOptVecSTEP(i/2+1) distOptVecSTEP(i/2+1)];
            end
            
            for i=2:2:2*N-1
                distOptACCEL(i)   = distOptACCEL(i-1) + d_f(i/2);
                distOptACCEL(i+1) = distOptACCEL(i) + d_ramps(i/2+1);
                vIterationACCEL = xOpt_accel(i/2);
                timeOptACCEL(i) = timeOptACCEL(i-1) + d_f(i/2)/vIterationACCEL;
                timeOptACCEL(i+1) = timeOptACCEL(i) + t_ramps(i/2);
            end
            
            distOptACCEL(end) = d_f(end) + distOptACCEL(end-1);
            timeOptACCEL(end) = d_f(end)/xOpt_accel(end) + timeOptACCEL(end-1);
            
end




