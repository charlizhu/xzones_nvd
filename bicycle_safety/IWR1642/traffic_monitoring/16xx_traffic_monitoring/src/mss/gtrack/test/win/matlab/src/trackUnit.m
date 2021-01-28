classdef trackUnit < matlab.System

    
properties (Nontunable)
    %trackUnitConfig, received from upper layer
    trackUnitConfig
end

properties (Nontunable)
    
    %gatingConstantVolume
    gatingConstantVolume = 16
    %gatingMaxLimits Absolute maximum limits: [range, width, velocity (default: 12m, 4m, nolimit]
    gatingMaxLimitRange = 8
    gatingMaxLimitWidth = 4
    gatingMaxLimitVelocity = 0

    %scoringDispersionFactor
    scoringDispersionFactor = 1
    
    %state changes thresholds, in frames
    detect2active = 3
    detect2free = 5
    active2free = 5
    static2free = 20
    exit2free = 3
    
    %Tracker default standard deviations for the measurment vector. Based on uniformly distributed RV in range [a b]: sqrt(1/12)*(b-a)
    % 1m range, 1m in width,  1m/s in doppler
    variationsRange = 1/sqrt(12)
    variationsWidth = 1/sqrt(12)
    variationsDoppler = 1/sqrt(12)
    
    %unrollingRangeRateAlpha
    unrollingRangeRateAlpha = 0.5
    %unrollingRangeRateConf, in %
    unrollingRangeRateConf = 0.1
end

properties (Access = private, Nontunable)
    pTid
    pHeartBeatCount
    pComputeFormat
    pStateVectorLength
    pMaxRadialVelocity
    pRadialVelocityResolution    
    pStateVectorType
    pMeasurementVectorLength
    pMaxNumPoints
    pMaxAcceleration
    pInitialRadialVelocity
    pMeasurementStd
    pColors

    pScenery
    
    pVariationsRange
    pVariationsWidth
    pVariationsDoppler
    
    pGatingConstantVolume
    pGatingMaxLimits
    
    pScoringDispersionFactor
    
    pVerbose
    pDeltaT

    state
    marker
    modelSize
    velocityHandling
    processVariance
    isTargetStatic
    allocationTime
    allocationRange
    allocationVelocity
    rangeRate
    tick
    color
    uHistory
    sHatHistory
    hSHistory
    ecHistory
    associatedPoints
    detect2activeCount
    detect2freeCount
    active2freeCount
    
    pDetect2active
    pDetect2free
    pActive2free
    pStatic2free
    pExit2free
    
    vehicleId

    Q
    F
    S_hat
    S_apriori_hat
    P
    P_apriori
    H_s_apriori_hat
    gD
    gC
    gC_inv
    gain
end

properties (Constant = true)
    ACTION_INSTANTIATE = 0;
    ACTION_PREDICT = 1;
    ACTION_UPDATE = 2;
end

properties (Constant = true)
    STATE_FREE = 0;
    STATE_INIT = 1;
    STATE_DETECTION = 2;
    STATE_ACTIVE = 3;
end


properties (Constant = true)
    FORMAT_2D = 0;
    FORMAT_2DA = 1;
    FORMAT_3D = 2;
end

properties (Constant = true)
    VELOCITY_INIT = 1;
    VELOCITY_RATE_FILTER = 2;
    VELOCITY_TRACKING = 3;
    VELOCITY_LOCKED = 4;
end

% veboseness masks
properties (Constant = true)
    VERBOSE_ERROR_INFO = 1;     %level 1
    VERBOSE_WARNING_INFO = 2;   %level 2
    VERBOSE_DEBUG_INFO = 4;     %level 3
    VERBOSE_MATRIX_INFO = 8;
    VERBOSE_UNROLL_INFO = 16;
    VERBOSE_STATE_INFO = 32;
    VERBOSE_ASSOSIATION_INFO = 64;
    VERBOSE_GATEXY_INFO = 128;
    VERBOSE_GATERA_INFO = 256;
    VERBOSE_GATEG1_INFO = 512;
end

methods
    
    function obj = trackUnit(varargin)
        setProperties(obj, nargin, varargin{:});
        
        obj.pScenery = obj.trackUnitConfig.advParams.scenery;
        
        if(~isempty(obj.trackUnitConfig.advParams.gating))
            obj.pGatingConstantVolume = obj.trackUnitConfig.advParams.gating(1);
            obj.pGatingMaxLimits = obj.trackUnitConfig.advParams.gating(2:4);
        else
            obj.pGatingConstantVolume = obj.gatingConstantVolume;
            obj.pGatingMaxLimits = [obj.gatingMaxLimitRange, obj.gatingMaxLimitWidth, obj.gatingMaxLimitVelocity];
        end
                
        if(~isempty(obj.trackUnitConfig.advParams.thresholds))
            obj.pDetect2active = obj.trackUnitConfig.advParams.thresholds(1);
            obj.pDetect2free = obj.trackUnitConfig.advParams.thresholds(2);
            obj.pActive2free = obj.trackUnitConfig.advParams.thresholds(3);
            obj.pStatic2free = obj.trackUnitConfig.advParams.thresholds(4);
            obj.pExit2free = obj.trackUnitConfig.advParams.thresholds(5);
        else         
            obj.pDetect2active = obj.detect2active;
            obj.pDetect2free = obj.detect2free;
            obj.pActive2free = obj.active2free;
            obj.pStatic2free = obj.static2free;
            obj.pExit2free = obj.exit2free;
        end
        
        if(~isempty(obj.trackUnitConfig.advParams.variations))
            obj.pVariationsRange = obj.trackUnitConfig.advParams.variations(1);
            obj.pVariationsWidth = obj.trackUnitConfig.advParams.variations(2);
            obj.pVariationsDoppler = obj.trackUnitConfig.advParams.variations(3);
        else         
            obj.pVariationsRange = obj.variationsHeight;
            obj.pVariationsWidth = obj.variationsWidth;
            obj.pVariationsDoppler = obj.variationsDoppler;
        end
        
        obj.pScoringDispersionFactor = obj.scoringDispersionFactor;

        obj.pMaxNumPoints = obj.trackUnitConfig.maxNumPoints;
        obj.pStateVectorType = obj.trackUnitConfig.stateVectorType;
        
        obj.pMaxRadialVelocity = obj.trackUnitConfig.maxRadialVelocity;
        obj.pRadialVelocityResolution = obj.trackUnitConfig.radialVelocityResolution;
        obj.pMaxAcceleration = obj.trackUnitConfig.maxAcceleration;
        obj.pInitialRadialVelocity = obj.trackUnitConfig.initialRadialVelocity;
        obj.pVerbose = obj.trackUnitConfig.verbose;
         
        obj.pDeltaT = obj.trackUnitConfig.deltaT;
        obj.pTid = obj.trackUnitConfig.tid;
        obj.pColors = 'bgrcm';
        
        dt = obj.pDeltaT;

        switch obj.pStateVectorType
            case obj.FORMAT_2D
                % We track 4 variables per track: S = [posx, posy, velx, vely]
                obj.pStateVectorLength = 4;
                % We expect 3 measurments: U = [range, azimuth, doppler]
                obj.pMeasurementVectorLength = 3;
                obj.pComputeFormat = obj.FORMAT_2D;
                
                obj.F = [
                    1 0 dt 0;
                    0 1 0 dt;
                    0 0 1 0;
                    0 0 0 1];
                
                a = 0.0001; b = 0.001; c = 0.01; d = 1;
                obj.Q = [a 0 b 0;
                         0 a 0 b;
                         b 0 c 0;
                         0 b 0 c]* d;                
                
            case obj.FORMAT_2DA
                % We track 6 variables per track: S = [posx posy velx vely accx accy]            
                obj.pStateVectorLength = 6;
                % We expect 3 measurments: U = [range, azimuth, doppler]
                obj.pMeasurementVectorLength = 3;
                obj.pComputeFormat = obj.FORMAT_2DA;
                obj.F = [
                    1 0 dt 0  dt^2/2 0;
                    0 1 0  dt 0      dt^2/2;
                    0 0 1  0  dt     0;
                    0 0 0  1  0      dt;
                    0 0 0  0  1      0;
                    0 0 0  0  0      1;
                    ];
                
                obj.Q = [
                    dt^4/4  0       dt^3/2  0       dt^2/2  0;
                    0       dt^4/4  0       dt^3/2  0       dt^2/2;
                    dt^3/2  0       dt^2    0       dt      0;
                    0       dt^3/2  0       dt^2    0       dt;
                    dt^2/2  0       dt      0       1       0;
                    0       dt^2/2  0       dt      0       1;
                    ];
                    
            otherwise
                error('radarTracking:invalidParameter', 'invalidStateTrackingType');
        end    
        obj.S_hat = zeros(obj.pStateVectorLength, 1);
        obj.S_apriori_hat = zeros(obj.pStateVectorLength,1);
        obj.H_s_apriori_hat = zeros(obj.pMeasurementVectorLength,1);
        obj.P = zeros(obj.pStateVectorLength, obj.pStateVectorLength);
        obj.P_apriori = zeros(obj.pStateVectorLength, obj.pStateVectorLength);
        obj.gD = zeros(obj.pMeasurementVectorLength, obj.pMeasurementVectorLength);
        obj.gC = zeros(obj.pMeasurementVectorLength, obj.pMeasurementVectorLength);
        obj.gC_inv = zeros(obj.pMeasurementVectorLength, obj.pMeasurementVectorLength);
%{        
        obj.state = obj.STATE_DETECTION;
%        obj.velocityHandling = obj.VELOCITY_INIT;
%        obj.modelSize = 4; %We always start with x,y,vx,vy
%        obj.processVariance = obj.initialProcessVariance;
        obj.velocityHandling = obj.VELOCITY_LOCKED;
        obj.modelSize = 6;
%}
        obj.allocationTime = 0;
        obj.vehicleId = 0;
        obj.rangeRate = 0;
        obj.marker = 0;
        obj.tick = 0;
        obj.gain = 0;
        obj.color = obj.pColors(mod(obj.pTid,length(obj.pColors))+1);
        obj.uHistory = zeros(400,3);
        obj.sHatHistory = zeros(400,6);
        obj.hSHistory = zeros(400,3);
        obj.ecHistory = zeros(400,3);
        obj.associatedPoints = 0;
        obj.detect2activeCount = 0;
        obj.detect2freeCount = 0;
        obj.active2freeCount = 0;
    end
    
    function trackUnitInitialize(obj, timeStamp, um)

        tid = obj.pTid;
        obj.pHeartBeatCount = timeStamp;
        obj.allocationTime = timeStamp;
        obj.allocationRange = um(1);
        obj.allocationVelocity = um(3);
        
        obj.state = obj.STATE_DETECTION;
        obj.isTargetStatic = false;
        obj.modelSize = 6;
        obj.processVariance = (1/2*obj.pMaxAcceleration).^2;
        obj.Q(1:2:end,:) = obj.Q(1:2:end,:) * obj.processVariance(1);
        obj.Q(2:2:end,:) = obj.Q(2:2:end,:) * obj.processVariance(2);
        
        obj.velocityHandling = obj.VELOCITY_INIT;
        rvIn = um(3);
        obj.rangeRate = unrollRadialVelocity(obj.pMaxRadialVelocity, obj.pInitialRadialVelocity, rvIn);

        
        % State Initialization
        um(3) = obj.rangeRate;
        obj.S_apriori_hat = computeCartesian(obj.pComputeFormat, um(1:3));
        obj.P_apriori = diag([0,0,0.5,0.5,1,1]);

        % For test purposes only
        obj.S_hat = computeCartesian(obj.pComputeFormat, um(1:3));
        obj.P = eye(6);

        obj.H_s_apriori_hat = um(1:3);
        if(bitand(obj.pVerbose, obj.VERBOSE_STATE_INFO))
            disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: ALLOCATED, Range ', num2str(um(1),3),', Angle ', num2str(um(2),3), ', Doppler ',num2str(rvIn,3), '=>', num2str(um(3),3)]);
        end
    end    
    
    function trackUnitPredict(obj)
        obj.pHeartBeatCount = obj.pHeartBeatCount + 1;
        % Prediction step of Kalman Filter
        
        tid = obj.pTid;
        mSize = obj.modelSize;

        if(obj.isTargetStatic)
            obj.S_apriori_hat(1:mSize) = obj.S_hat(1:mSize);
            obj.P_apriori(1:mSize,1:mSize) = obj.P(1:mSize,1:mSize);
        else
            obj.S_apriori_hat(1:mSize) = obj.F(1:mSize,1:mSize) * obj.S_hat(1:mSize);
            obj.P_apriori(1:mSize,1:mSize) = obj.F(1:mSize,1:mSize) * obj.P(1:mSize,1:mSize) * obj.F(1:mSize,1:mSize)' + obj.Q(1:mSize,1:mSize);
            % Enforce symmetry constraint on P_apriori
            obj.P_apriori(1:mSize,1:mSize) = 1/2 * (obj.P_apriori(1:mSize,1:mSize)+obj.P_apriori(1:mSize,1:mSize)');
        end
            
        % Convert from spherical to Cartesian coordinates
        obj.H_s_apriori_hat = computeH(obj.pComputeFormat, obj.S_apriori_hat);
        
        if(bitand(obj.pVerbose,obj.VERBOSE_MATRIX_INFO))
            disp(['[',num2str(tid),']: Predict: s-hat=>s-apriori']); disp([obj.S_hat,obj.S_apriori_hat]);
            disp(['[',num2str(tid),']: Predict: P,Q => P-apriori']); disp(obj.P(1:mSize,1:mSize)); disp(obj.Q(1:mSize,1:mSize)), disp(obj.P_apriori(1:mSize,1:mSize));
            disp(['[',num2str(tid),']: Predict: H-s =']); disp(obj.H_s_apriori_hat);
        end
    end
    
    function [scoreOut, indexOut, rvOut, centroid, ec, G] = trackUnitAssociate(obj, u, scoreIn, indexIn)
        urv = u(1:3,:);
        nMeas = size(u,2);
        rvOut = zeros(nMeas,1);
        scoreOut = scoreIn;
        indexOut = indexIn;

        % Compute ellipsoid expansion factor G
        if(obj.isTargetStatic)
            G = 1; a = 0; b = 0; c = 0;
        else
            [G, a, b, c] = gateCreateLim(obj.pGatingConstantVolume, obj.gC_inv, obj.H_s_apriori_hat, obj.pGatingMaxLimits);
        end
        if(bitand(obj.pVerbose, obj.VERBOSE_MATRIX_INFO))
            disp(['[',num2str(obj.pTid),']: G,[a,b,c]']); disp([G,a,b,c]);
        end
        a = det(obj.gC);
        hs = obj.H_s_apriori_hat;

        if(obj.velocityHandling < obj.VELOCITY_LOCKED)
            hs(3) = obj.allocationVelocity;
        end
        
        urv(3,:) = unrollRadialVelocity(obj.pMaxRadialVelocity, hs(3), u(3,:));
        
        for k = 1:nMeas
            
            % Gating step
            % Compute the amount of innovation this measurement
            % introduces wrt the tracking object
            u_tilda = urv(:,k) - hs;
            chi2 = u_tilda'*obj.gC_inv*u_tilda;
            if(chi2 < G)
                % Compute scoring function
                score = obj.pScoringDispersionFactor*log(a)+chi2;
                if(score < scoreIn(k))
                    scoreOut(k) = score;
                    indexOut(k) = obj.pTid;
                    rvOut(k) = urv(3,k);
                end
            end
        end
        centroid = obj.H_s_apriori_hat;
        ec = obj.gC_inv;
        obj.gain = G;
    end
    
    function state = trackUnitUpdate(obj, u, ~, mid)

        tid = obj.pTid;
        % Measurement update
        % set the model size
        mSize = obj.modelSize;
        
        num = length(mid);
      
        obj.associatedPoints = obj.associatedPoints + num;
        
        if(num == 0)
           % No measurments available
           % Set state and process based on prediction
           obj.S_hat = obj.S_apriori_hat;
           obj.P = obj.P_apriori;
           % Check whether target is static
           if((abs(obj.S_hat(3)) < 3*obj.pRadialVelocityResolution) && (abs(obj.S_hat(4)) < 3*obj.pRadialVelocityResolution))
               % Target is static => force zero velocity and zero acceleration           
               obj.S_hat(3:end) = 0;
               % Zero down process variance
               obj.isTargetStatic = true;
           end
           % Report a miss event
           state = trackEvent(obj, 0, 0);
           return;
           
        else
            obj.isTargetStatic = false;
            % Getting back to dynamic model 
        end
        
        um = mean(u(1:2,mid),2);
        
        velPilot = u(3,mid(1));
        u(3,mid(2:end)) = unrollRadialVelocity(obj.pMaxRadialVelocity, velPilot, u(3,mid(2:end)));
        velMeasured = mean(u(3,mid));
   
        switch(obj.velocityHandling)

            case obj.VELOCITY_INIT
                um(3) = obj.rangeRate;
                obj.velocityHandling = obj.VELOCITY_RATE_FILTER;

                if(bitand(obj.pVerbose, obj.VERBOSE_UNROLL_INFO))
                    disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update State VINIT=>VFILT, ', num2str(velMeasured,3), '=>', num2str(um(3),3)]);
                end


            case obj.VELOCITY_RATE_FILTER
                instanteneousRangeRate = (um(1)-obj.allocationRange)/((obj.pHeartBeatCount-obj.allocationTime)*obj.pDeltaT);
                obj.rangeRate = obj.unrollingRangeRateAlpha*obj.rangeRate + (1-obj.unrollingRangeRateAlpha)*instanteneousRangeRate;

                um(3) = unrollRadialVelocity(obj.pMaxRadialVelocity, obj.rangeRate, velMeasured);

                rrError = (instanteneousRangeRate-obj.rangeRate)/obj.rangeRate;
                if(abs(rrError) < obj.unrollingRangeRateConf)
                    % RangeRate is stable enough to use it against H-s
                    if(bitand(obj.pVerbose,obj.VERBOSE_UNROLL_INFO))
                        disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update State VRFIL => VTRACK, Unrolling with RangeRate ', num2str(obj.rangeRate,3),', H-s ', num2str(obj.H_s_apriori_hat(3),3), ', ', num2str(velMeasured,3),'=>',num2str(um(3),3)]);
                    end
                    obj.velocityHandling = obj.VELOCITY_TRACKING;
                else
                    % RangeRate is stable enough to use it against H-s
                    if(bitand(obj.pVerbose,obj.VERBOSE_UNROLL_INFO))
                        disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update State VRFIL, Unrolling with RangeRate ', num2str(obj.rangeRate,3),', H-s ', num2str(obj.H_s_apriori_hat(3),3), ', ', num2str(velMeasured,3),'=>',num2str(um(3),3)]);
                    end
                end

            case obj.VELOCITY_TRACKING
                instanteneousRangeRate = (um(1)-obj.allocationRange)/((obj.pHeartBeatCount-obj.allocationTime)*obj.pDeltaT);
                obj.rangeRate = obj.unrollingRangeRateAlpha*obj.rangeRate + (1-obj.unrollingRangeRateAlpha)*instanteneousRangeRate;

                um(3) = unrollRadialVelocity(obj.pMaxRadialVelocity, obj.rangeRate, velMeasured);

                velError = abs(obj.H_s_apriori_hat(3) - um(3));
                if (velError < abs(um(3))/10)
                    if(bitand(obj.pVerbose,obj.VERBOSE_UNROLL_INFO))
                        disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update State VTRACK => VLOCK, Unrolling with RangeRate ', num2str(obj.rangeRate,3),', H-s ', num2str(obj.H_s_apriori_hat(3),3), ', ', num2str(velMeasured,3),'=>',num2str(um(3),3)]);
                    end
                    obj.velocityHandling = obj.VELOCITY_LOCKED;
                else
                    if(bitand(obj.pVerbose,obj.VERBOSE_UNROLL_INFO))
                        disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update State VTRACK, Unrolling with RangeRate ', num2str(obj.rangeRate,3),', H-s ', num2str(obj.H_s_apriori_hat(3),3), ', ', num2str(velMeasured,3),'=>',num2str(um(3),3)]);
                    end
                end
                
            case obj.VELOCITY_LOCKED
                um(3) = unrollRadialVelocity(obj.pMaxRadialVelocity, obj.H_s_apriori_hat(3), velMeasured);
        end
        
        % Compute measurement noise covariance matrix, Rg
        rangeVar = obj.pVariationsRange*obj.pVariationsRange;
        angleVar = (2*atan(0.5*obj.pVariationsWidth/obj.H_s_apriori_hat(1)))^2;
        doplperVar = obj.pVariationsDoppler*obj.pVariationsDoppler;

        Rm = diag([rangeVar angleVar doplperVar]);
        D = zeros(mSize,mSize);
        
        % Update group dispersion matrix
        if(num > 3)
            % New dispersion is computed only for large sets
            d11 = mean((u(1,mid)-um(1)).^2);
            d22 = mean((u(2,mid)-um(2)).^2);
            d33 = mean((u(3,mid)-velMeasured).^2);

            d12 = mean((u(1,mid)-um(1)).*(u(2,mid)-um(2)));
            d13 = mean((u(1,mid)-um(1)).*(u(3,mid)-velMeasured));
            d23 = mean((u(2,mid)-um(2)).*(u(3,mid)-velMeasured));
            
            D = [d11,d12,d13; d12,d22,d23; d13,d23,d33];
            
            if(obj.associatedPoints > 50)
                alpha = num/50;
            else
                alpha = num/obj.associatedPoints;
            end
            % Update
            obj.gD = (1-alpha)*obj.gD + alpha*D;
        end
        
        % Compute Measurement noise covariance matrix Rc, used for Kalman
        % updates
        if(num > 10)
            alpha = 0;
        else
            alpha = (10-num)/((10-1)*num);
        end
            
        Rc = Rm/num + alpha*obj.gD*eye(3);
        
        % Compute the Jacobian (can be over-written)
        J = computeJacobian(obj.pComputeFormat, obj.S_apriori_hat);

        % a) innovation (or residual vector in [1])
        u_tilda = um - obj.H_s_apriori_hat;

        % b) innovation covariance (or residual covariance) of the centroid
        cC = J(:,1:mSize) * obj.P_apriori(1:mSize,1:mSize) * J(:,1:mSize)' + Rc;
        cC_inv = matinv(cC);

        obj.gC =  J(:,1:mSize) * obj.P_apriori(1:mSize,1:mSize) * J(:,1:mSize)' + Rm + obj.gD;
        obj.gC_inv = matinv(obj.gC);

        if(bitand(obj.pVerbose, obj.VERBOSE_DEBUG_INFO))
            disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update U={', num2str(um(1),3),', ',num2str(um(2),3),', ', num2str(velMeasured,3),'=>',num2str(um(3),3),'}']);
        end

        % c) Kalman gain using matrix inverse via Cholesky decomposition
        K = obj.P_apriori(1:mSize,1:mSize) * J(:,1:mSize)' * cC_inv;
        % d) State estimation
        obj.S_hat(1:mSize) = obj.S_apriori_hat(1:mSize) + K * u_tilda;
        % e) Estimation covariance
        obj.P(1:mSize,1:mSize) = obj.P_apriori(1:mSize,1:mSize) - K * J(:,1:mSize) * obj.P_apriori(1:mSize,1:mSize);
        
        if(bitand(obj.pVerbose, obj.VERBOSE_MATRIX_INFO))
            disp('Rm = '); disp(Rm);
            disp('D = '); disp(D);
            disp('gD = '); disp(obj.gD);
            disp('Rc = '); disp(Rc);
            
            disp('J = '); disp(J(:,1:mSize));
            disp('u_tilda ='); disp(u_tilda');
            
            disp('P = '); disp(obj.P_apriori);
            disp('cC = '); disp(cC);
            disp('cC_inv = '); disp(cC_inv);

            disp('K ='); disp(K);            

            disp(['[',num2str(tid),']: Update: s-apriori=>s-hat']); disp([obj.S_apriori_hat, obj.S_hat]);
            disp('P_hat = '); disp(obj.P(1:mSize,1:mSize));

            disp('gC = '); disp(obj.gC);
            disp('gC_inv = '); disp(obj.gC_inv);

        end
        
        if(bitand(obj.pVerbose, obj.VERBOSE_DEBUG_INFO))
            disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update S={',num2str(obj.S_hat(1),3),', ',num2str(obj.S_hat(2),3),', ',...
                num2str(obj.S_hat(3),3),', ',num2str(obj.S_hat(4),3),', ',num2str(obj.S_hat(5),3),', ',num2str(obj.S_hat(6),3),'}']);
        end
        
        %update state machine with hit event
        state = trackEvent(obj, num, um);      
    end
    
    function [S, state, ec, g] = trackUnitGetInfo(obj)
        state = obj.state;
        S = obj.S_hat;
        ec = obj.gC_inv;
        g = obj.gain;
    end
    
    function isBehind = trackUnitIsBehind(obj, um)
        if( (abs(um(2) - obj.H_s_apriori_hat(2)) < 2*pi/180) && ...
                (um(1) > obj.H_s_apriori_hat(1)) )  
            isBehind = true;
        else
            isBehind = false;
        end
    end
    
    function trackUnitSetMarker(obj, value)
        obj.marker = value;
    end
    
    function value = trackUnitGetMarker(obj)
        value = obj.marker;
    end

end

methods (Access = protected)
       
    function setupImpl(obj)
        obj.pGatingConstantVolume = obj.gatingConstantVolume;
        obj.pGatingMaxLimits = obj.gatingMaxLimits;
        obj.pScoringDispersionFactor = obj.scoringDispersionFactor;
    end            

    function state = trackEvent(obj, num, um)
        
        tid = obj.pTid;
        
        % Collect tracker states history
        ind = mod(obj.tick, 400)+1; % one-based
        if(num)
            obj.uHistory(ind,:) = um;
        else
            obj.uHistory(ind,:) = zeros(1,3);
        end
        obj.sHatHistory(ind,:) = obj.S_hat;
        obj.hSHistory(ind,:) = obj.H_s_apriori_hat;
        obj.ecHistory(ind,:) = [obj.gC(1,1), obj.gC(2,2), obj.gC(3,3)];
        obj.tick = obj.tick + 1;
        
        switch obj.state
            case obj.STATE_DETECTION
                if(num > 5)
                    % Hit
                    obj.detect2freeCount = 0;
                    obj.detect2activeCount = obj.detect2activeCount +1;
                    if(obj.detect2activeCount > obj.pDetect2active)
                        if(bitand(obj.pVerbose, obj.VERBOSE_STATE_INFO))
                            disp(['[',num2str(tid),']: DET => ACTIVE']);
                        end
                        obj.state = obj.STATE_ACTIVE;
                    end
                else
                    if(num == 0)
                        % Miss
                        obj.detect2freeCount = obj.detect2freeCount +1;
                        if(obj.detect2activeCount > 0)
                            obj.detect2activeCount = obj.detect2activeCount -1;
                        end
                        
                        if(obj.detect2freeCount > obj.pDetect2free)
                            if(bitand(obj.pVerbose, obj.VERBOSE_STATE_INFO))
                                disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: DET => FREE']);
                            end
                            obj.state = obj.STATE_FREE;
                        end
                    end
                end
            case obj.STATE_ACTIVE
                if(num)
                    % Hit
                    obj.active2freeCount = 0;
                else
                    % Miss
                    obj.active2freeCount = obj.active2freeCount + 1;
                    inBoxInd = false;
                    for n = 1:obj.pScenery.numStaticBoxes
                        box = obj.pScenery.staticBoxes(n,:);
                        inBoxInd = inBoxInd || ((obj.S_hat(1) > box(1)) && (obj.S_hat(1) < box(2)) && (obj.S_hat(2) > box(3)) && (obj.S_hat(2) < box(4)));            
                    end
                    if ~inBoxInd
                        % If we are in entry/exit zone
                        thre = obj.pExit2free;
                    elseif (obj.isTargetStatic)
                        % If we are static    
                        thre = obj.pStatic2free; 
                    else
                        thre = obj.pActive2free;
                    end
                    if(thre > obj.tick)
                        thre = obj.tick;
                    end
                    
                    if(obj.active2freeCount > thre)
                        if(bitand(obj.pVerbose, obj.VERBOSE_STATE_INFO))
                            disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: ACTIVE => FREE']);
                        end
                        obj.state = obj.STATE_FREE;
                    end
                end
        end
        state = obj.state;
    end
end
end