classdef trackUnit < matlab.System

    
properties (Nontunable)
    %trackUnitConfig, received from upper layer
    trackUnitConfig
end

properties (Nontunable)
    %gatingConstantVolume
    gatingConstantVolume = 16
    %gatingMaxLimits Absolute maximum limits: [range, angle, velocity]
    gatingMaxLimits = [12, 30*pi/180, 0]
    %scoringDispersionFactor
    scoringDispersionFactor = 1
    %trackerActiveThreshold in frames
    trackerActiveThreshold = 3
    %trackerForgetThreshold
    trackerForgetThreshold = 3
    %allocationThresholdPoints, minimum number of points
    allocationThresholdPoints = 3    
    %allocationThresholdDistance, in meters
    allocationThresholdDistance = 4
    %allocationThresholdVelocity, in m/c
    allocationThresholdVelocity = 2
    %unrollingRangeRateAlpha
    unrollingRangeRateAlpha = 0.5
    %unrollingRangeRateConf, in %
    unrollingRangeRateConf = 0.1
    %initialProcessVariance
    initialProcessVariance = 1000
end

properties (Access = private, Nontunable)
    pTid
    pHeartBeatCount
    pComputeFormat
    pStateVectorLength
    pMaxURadialVelocity
    pStateVectorType
    pMeasurementVectorLength
    pMaxNumPoints
    pMaxAcceleration
    pColors
    
    pGatingConstantVolume
    pGatingMaxLimits
    pScoringDispersionFactor
    pTrackerActiveThreshold
    pTrackerForgetThreshold
    pVerbose
    pDeltaT

    state
    marker
    modelSize
    velocityHandling
    processVariance
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
    activeThreshold
    forgetThreshold
    vehicleId

    Q
    F
    S_hat
    S_apriori_hat
    P
    P_apriori
    H_s_apriori_hat
    iCgroup_inv
    iCgroup
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
        
        obj.pGatingConstantVolume = obj.gatingConstantVolume;
        obj.pGatingMaxLimits = obj.gatingMaxLimits;
        obj.pScoringDispersionFactor = obj.scoringDispersionFactor;

        obj.pMaxNumPoints = obj.trackUnitConfig.maxNumPoints;
        obj.pStateVectorType = obj.trackUnitConfig.stateVectorType;
        
        obj.pMaxURadialVelocity = obj.trackUnitConfig.maxURadialVelocity;
        obj.pMaxAcceleration = obj.trackUnitConfig.maxAcceleration;
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
        obj.iCgroup = zeros(obj.pMeasurementVectorLength, obj.pMeasurementVectorLength);
        obj.iCgroup_inv = zeros(obj.pMeasurementVectorLength, obj.pMeasurementVectorLength);
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
        obj.color = obj.pColors(mod(obj.pTid,length(obj.pColors))+1);
        obj.uHistory = zeros(400,3);
        obj.sHatHistory = zeros(400,6);
        obj.hSHistory = zeros(400,3);
        obj.ecHistory = zeros(400,3);
        obj.associatedPoints = 1;
        obj.detect2activeCount = 0;
        obj.detect2freeCount = 0;
        obj.active2freeCount = 0;
        obj.activeThreshold = obj.trackerActiveThreshold;
        obj.forgetThreshold = obj.trackerForgetThreshold;
    end
    
    function trackUnitInitialize(obj, timeStamp, um, vid)

        tid = obj.pTid;
        obj.pHeartBeatCount = timeStamp;
        obj.allocationTime = timeStamp;
        obj.allocationRange = um(1);
        obj.allocationVelocity = um(3);
        obj.vehicleId = vid;
        
        obj.state = obj.STATE_DETECTION;
        obj.modelSize = 6;
        obj.processVariance = (1/2*obj.pMaxAcceleration)^2;
        
        obj.velocityHandling = obj.VELOCITY_INIT;
        rvIn = um(3);
        obj.rangeRate = unrollRadialVelocity(obj.pMaxURadialVelocity,-20, rvIn);

        
        % State Initialization
        um(3) = obj.rangeRate;
        obj.S_apriori_hat = computeCartesian(obj.pComputeFormat, um, 0);
        obj.P_apriori = eye(6);

        % For test purposes only
        obj.S_hat = computeCartesian(obj.pComputeFormat, um, 0);
        obj.P = eye(6);

        obj.H_s_apriori_hat = um;
        if(bitand(obj.pVerbose, obj.VERBOSE_STATE_INFO))
            disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: ALLOCATED, Range ', num2str(um(1),3),', Angle ', num2str(um(2),3), ', Doppler ',num2str(rvIn,3), '=>', num2str(um(3),3)]);
        end
    end    
    
    function trackUnitPredict(obj)
        obj.pHeartBeatCount = obj.pHeartBeatCount + 1;
        % Prediction step of Kalman Filter
        
        tid = obj.pTid;
        mSize = obj.modelSize;
        
        obj.S_apriori_hat(1:mSize) = obj.F(1:mSize,1:mSize) * obj.S_hat(1:mSize);
        obj.P_apriori(1:mSize,1:mSize) = obj.F(1:mSize,1:mSize) * obj.P(1:mSize,1:mSize) * obj.F(1:mSize,1:mSize)' + obj.Q(1:mSize,1:mSize)*obj.processVariance;
        % Enforce symmetry constraint on P_apriori
        obj.P_apriori(1:mSize,1:mSize) = 1/2 * (obj.P_apriori(1:mSize,1:mSize)+obj.P_apriori(1:mSize,1:mSize)');
        % Convert from spherical to Cartesian coordinates
        obj.H_s_apriori_hat = computeH(obj.pComputeFormat, obj.S_apriori_hat);
        
        if(bitand(obj.pVerbose,obj.VERBOSE_MATRIX_INFO))
            disp(['[',num2str(tid),']: Predict: s-hat=>s-apriori']); disp([obj.S_hat,obj.S_apriori_hat]);
            disp(['[',num2str(tid),']: Predict: P,Q => P-apriori']); disp(obj.P(1:mSize,1:mSize)); disp(obj.Q(1:mSize,1:mSize)*obj.processVariance), disp(obj.P_apriori(1:mSize,1:mSize));
            disp(['[',num2str(tid),']: Predict: H-s =']); disp(obj.H_s_apriori_hat);
        end
    end
    
    function [scoreOut, indexOut, rvOut, centroid, ec, G] = trackUnitAssociate(obj, u, scoreIn, indexIn)
        urv = u;
        nMeas = size(u,2);
        tHistStruct = struct('num', 0, 'rvExp', 0, 'rvSum', 0, 'mid', [], 'rvOut', [], 'score', []);
        bucket = repmat(tHistStruct, 1, 4);
        rvOut = zeros(nMeas,1);
        scoreOut = scoreIn;
        indexOut = indexIn;
        
        % Compute ellipsoid expansion factor G
        [G, a, b, c] = gateCreateLim(obj.pGatingConstantVolume, obj.iCgroup_inv, obj.pGatingMaxLimits);        
        if(bitand(obj.pVerbose, obj.VERBOSE_MATRIX_INFO))
            disp(['[',num2str(obj.pTid),']: G,[a,b,c]']); disp([G,a,b,c]);
        end
        a = det(obj.iCgroup);
        hs = obj.H_s_apriori_hat;

        if(obj.velocityHandling < obj.VELOCITY_TRACKING)
            hs(3) = obj.allocationVelocity;
        end
        
        urv(3,:) = unrollRadialVelocity(obj.pMaxURadialVelocity, hs(3), u(3,:));
        
        for k = 1:nMeas
            % Gating step
            % Compute the amount of innovation this measurement
            % introduces wrt the tracking object
            u_tilda = urv(:,k) - hs;
            chi2 = u_tilda'*obj.iCgroup_inv*u_tilda;
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
%{            
            
            else
            % Don't have stable radial velocity estimation yet
            rvExpKnown = 0;
            maxBucketCount = 1;
            maxBucketIndex = 1;
            for k = 1:nMeas
                % Gating step with no velocity data
                rvIn = u(3,k);
                u_tilda = u(:,k) - obj.H_s_apriori_hat;
                u_tilda(3) = 0;
                chi2 = u_tilda'*obj.iCgroup_inv*u_tilda;
                if(chi2 < G)
                    % Compute scoring function
                    score = obj.pScoringDispersionFactor*log(a)+chi2;

                    if(rvExpKnown == 1)
                        rv = unrollRadialVelocity(obj.pMaxURadialVelocity, rvExp, rvIn);
                        if(abs(rv - rvExp) < 2)
                            if(score < scoreIn(k))
                                scoreOut(k) = score;
                                indexOut(k) = obj.pTid;
                                rvOut(k) = rv;
                            end
                        end
                    else
                        % Place into one of four bins
                        for b=1:4
                            if(bucket(b).num == 0)
                                bucket(b).rvExp = rvIn;
                                bucket(b).rvSum = rvIn;
                                bucket(b).num = 1;
                                bucket(b).mid = k;
                                bucket(b).rvOut = rvIn;
                                bucket(b).score = score;
                                break;
                            else
                                rv = unrollRadialVelocity(obj.pMaxURadialVelocity, bucket(b).rvExp, rvIn);
                                if(abs(rv - bucket(b).rvExp) < 2)
                                    % Found a bucket
                                    bucket(b).rvSum = bucket(b).rvSum + rv;
                                    num = bucket(b).num + 1;
                                    bucket(b).rvExp = bucket(b).rvSum/num;
                                    bucket(b).mid = [bucket(b).mid k];
                                    bucket(b).rvOut = [bucket(b).rvOut rv];
                                    bucket(b).score = [bucket(b).score score];
                                    bucket(b).num = num;
                                    if(num > maxBucketCount)
                                        maxBucketCount = num;
                                        maxBucketIndex = b;
                                        if(num > 4)
                                            rvExpKnown = 1;
                                            rvExp = bucket(b).rvExp;
                                        end
                                    end
                                    break;
                                end
                            end
                        end                                
                    end
                end
            end
            
            for b = 1:bucket(maxBucketIndex).num
                score = bucket(maxBucketIndex).score(b);
                k = bucket(maxBucketIndex).mid(b);
                if(score < scoreIn(k))   
                    scoreOut(k) = score;
                    indexOut(k) = obj.pTid;            
                    rvOut(k) = bucket(maxBucketIndex).rvOut(b);         
                end
            end
        end
%}

        centroid = obj.H_s_apriori_hat;
        ec = obj.iCgroup_inv;
    end
    
    function state = trackUnitUpdate(obj, u, uvar, mid)

        tid = obj.pTid;
        num = length(mid);
        if(num == 0)
           % No measurments available, set estimations based on predictions */
           obj.S_hat = obj.S_apriori_hat;
           obj.P = obj.P_apriori;
           
           state = trackEvent(obj, 0, 0);
           return;
        end
        um = mean(u(1:2,mid),2);
        
        velPilot = u(3,mid(1));
        u(3,mid(2:end)) = unrollRadialVelocity(obj.pMaxURadialVelocity, velPilot, u(3,mid(2:end)));
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

                um(3) = unrollRadialVelocity(obj.pMaxURadialVelocity, obj.rangeRate, velMeasured);

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
                % Unroll radial veloctiy with H_s
                um(3) = unrollRadialVelocity(obj.pMaxURadialVelocity, obj.H_s_apriori_hat(3), velMeasured);

                velError = abs(obj.H_s_apriori_hat(3) - um(3));
                if (velError < abs(um(3))/10)
                    if(bitand(obj.pVerbose,obj.VERBOSE_UNROLL_INFO))
                        disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update State VTRACK => VLOCK, Unrolling with H-s ', num2str(obj.H_s_apriori_hat(3),3), ', ', num2str(velMeasured,3),'=>',num2str(um(3),3)]);
                    end
                    obj.velocityHandling = obj.VELOCITY_LOCKED;
                else
                    if(bitand(obj.pVerbose,obj.VERBOSE_UNROLL_INFO))
                        disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update State VTRACK, Unrolling with H-s = ', num2str(obj.H_s_apriori_hat(3),3), ': ', num2str(velMeasured,3),'=>',num2str(um(3),3)]);
                    end
                end
                
            case obj.VELOCITY_LOCKED
                um(3) = unrollRadialVelocity(obj.pMaxURadialVelocity, obj.H_s_apriori_hat(3), velMeasured);
        end
        
        % Calculate centroid variances
        uvarm = sum(uvar(:,mid),2)./(num^2);
        % variance can not be less than quantization
        uvarm = max([(0.4/2)^2;(pi/180/2)^2;(0.3298/2)^2], uvarm);
        
        R = diag(uvarm);
        
        % Compute dispersion matrix S
        s11 = mean((u(1,mid)-um(1)).^2);
        s22 = mean((u(2,mid)-um(2)).^2);
        s33 = mean((u(3,mid)-velMeasured).^2);
        s12 = mean((u(1,mid)-um(1)).*(u(2,mid)-um(2)));
        s13 = mean((u(1,mid)-um(1)).*(u(3,mid)-velMeasured));
        s23 = mean((u(2,mid)-um(2)).*(u(3,mid)-velMeasured));
        S = [s11,s12,s13;s12,s22,s23; s13,s23,s33];
        
        % Velocity handling
        % um(3) = unrollRadialVelocity(obj.pMaxURadialVelocity,obj.H_s_apriori_hat(3), um(3));
        
        % Enforce minimums
        % S = max(S,[0.4*0.4, 0.4*2*pi/180, 0.4*0.3298; 2*pi/180*0.4, 2*pi/180*2*pi/180, 2*pi/180*0.3298; 0.3298*0.4, 0.3298*2*pi/180, 0.3298*0.3298]);
        % S = max(S,[0, 0, 0; 0, 0, 0; 0, 0, 0.3298*0.3298]);


        % Compute the Jacobian (can be over-written)
        J = computeJacobian(obj.pComputeFormat, obj.S_apriori_hat);

        % Measurement update
        % set the model size
        mSize = obj.modelSize;

        % a) innovation (or residual vector in [1])
        u_tilda = um - obj.H_s_apriori_hat;

        % b) innovation covariance (or residual covariance)
        iC = J(:,1:mSize) * obj.P_apriori(1:mSize,1:mSize) * J(:,1:mSize)' + R;
        iC_inv = matinv(iC);

        obj.iCgroup = iC + S;
        obj.iCgroup_inv = matinv(obj.iCgroup);

        if(bitand(obj.pVerbose, obj.VERBOSE_DEBUG_INFO))
            disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update U={', num2str(um(1),3),', ',num2str(um(2),3),', ', num2str(velMeasured,3),'=>',num2str(um(3),3),'}']);
        end

        % c) Kalman gain using matrix inverse via Cholesky decomposition
        K = obj.P_apriori(1:mSize,1:mSize) * J(:,1:mSize)' * iC_inv;
        % d) State estimation
        obj.S_hat(1:mSize) = obj.S_apriori_hat(1:mSize) + K * u_tilda;
        % e) Estimation covariance
        obj.P(1:mSize,1:mSize) = obj.P_apriori(1:mSize,1:mSize) - K * J(:,1:mSize) * obj.P_apriori(1:mSize,1:mSize);
        
        if(bitand(obj.pVerbose, obj.VERBOSE_MATRIX_INFO))
            disp('R = '); disp(R);
            disp('S = '); disp(S);
            disp('J = '); disp(J(:,1:mSize));
            disp('u_tilda ='); disp(u_tilda');
            
            disp('iC = '); disp(iC);
            disp('iC_inv = '); disp(iC_inv);

            disp('iCgroup = '); disp(obj.iCgroup);
            disp('iCgroup_inv = '); disp(obj.iCgroup);

            disp('K ='); disp(K);            
            disp(['[',num2str(tid),']: Update: s-apriori=>s-hat']); disp([obj.S_apriori_hat, obj.S_hat]);
            disp('P_hat = '); disp(obj.P(1:mSize,1:mSize));
        end
        
        if(bitand(obj.pVerbose, obj.VERBOSE_DEBUG_INFO))
            disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: Update S={',num2str(obj.S_hat(1),3),', ',num2str(obj.S_hat(2),3),', ',...
                num2str(obj.S_hat(3),3),', ',num2str(obj.S_hat(4),3),', ',num2str(obj.S_hat(5),3),', ',num2str(obj.S_hat(6),3),'}']);
        end
        
        %update state machine with hit event
        state = trackEvent(obj, num, um);      
    end
    
    function [S, state] = trackUnitGetInfo(obj)
        state = obj.state;
        S = obj.S_hat;
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
    
        obj.associatedPoints = num;
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
        obj.ecHistory(ind,:) = [obj.iCgroup(1,1), obj.iCgroup(2,2), obj.iCgroup(3,3)];
        obj.tick = obj.tick + 1;
        
        switch obj.state
            case obj.STATE_DETECTION
                if(num > 5)
                    % Hit
                    obj.detect2freeCount = 0;
                    obj.detect2activeCount = obj.detect2activeCount +1;
                    if(obj.detect2activeCount > obj.activeThreshold)
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
                        
                        if(obj.detect2freeCount > obj.forgetThreshold)
                            if(bitand(obj.pVerbose, obj.VERBOSE_STATE_INFO))
                                disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: DET => FREE']);
                            end
                            obj.state = obj.STATE_FREE;
                        end
                    end
                end
            case obj.STATE_ACTIVE
                if(num > 5)
                    % Hit
                    if(obj.active2freeCount > 0)
                        obj.active2freeCount = obj.active2freeCount -1;
                    end                
                else
                    if (num == 0)    
                        % Miss
                        obj.active2freeCount = obj.active2freeCount +1;
                        if(obj.active2freeCount > obj.forgetThreshold)
                            if(bitand(obj.pVerbose, obj.VERBOSE_STATE_INFO))
                                disp([num2str(obj.pHeartBeatCount),': TID[',num2str(tid-1),']: ACTIVE => FREE']);
                            end
                            obj.state = obj.STATE_FREE;
                        end
                    end
                end
        end
        state = obj.state;
    end
end
end