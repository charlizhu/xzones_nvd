classdef trackModule < matlab.System

%radarTracking is Radar Target Tracking module
%   H = radarTracking creates a trackig module H which is an instance of
%   the radarTracking class. Tracking module estimates the position and
%   velocity of multiple radar targets in noisy measurement environment. It
%   employs Kalman filter for tracking object position and velocity [2][1].
%   Module expects measurement results which include estimation of range,
%   azimuth angle, elevation angle, and target velocity for each of the
%   detected targets.
%
%   H = radarTracking(Name,Value) creates a tracking object, H, with the
%   specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Module uses Kalman filter to “refine” the location estimates. Module
%   supports 2D or 3D target modeling The state of the Kalman filter at
%   time instant n is 4x1 (in 2D) or 6x1 (in 3D) vector defined as:
%       s(n)= [x(n) y(n) vx(n) vy(n) ]', (2D) 
%       s(n)= [x(n) y(n) z(n) vx(n) vy(n) vz(n)]', (3D) 
%       where
%           x(n), y(n), and z(n) are target positions, and 
%           vx(n), vy(n), and vz(n) are velocity projections
%
%   State transition is modeled as:
%       s(n)=F*s(n-1)+w(n), 
%       where
%           F is a state transition 4x4 (2D) 6x6 (2D with Acceleration) or 6x6 (3D) matrix. We define
%           F = [1 0 dt 0;
%                0 1 0 dt;
%                0 0 1 0;
%                0 0 0 1]; (2D)
%           F = [1  0  dt 0  dt^2/2  0;
%                0  1  0  dt 0  dt^2/2;
%                0  0  1  0  dt 0;
%                0  0  0  1  0  dt;
%                0  0  0  0  1  0;
%                0  0  0  0  0  1]; (2DA), and
%           F = [1  0  0  dt 0  0;
%                0  1  0  0  dt 0;
%                0  0  1  0  0  dt;
%                0  0  0  1  0  0;
%                0  0  0  0  1  0;
%                0  0  0  0  0  1], (3D).
%           w(n) is the vector of process noise with covariance matrix Q(n)
%           of size 4x4 (2D) or 6x6 (2DA), and 6x6 (3D).
%
%   The measuremnt vector u(n) provides Kalman filter with updates on
%   target location and velocity is defined as:
%       u(n) = [range(n) azimuth(n) velocity(n)] (2D)
%       u(n) = [range(n) azimuth(n) elevation(n) velocity(n)] (3D).
%
%   The relation between measurements updates and the state of the Kalman
%   filter is non-linear, and is defined as:
%       u(n) = H(S(n)) + v(n), 
%       where
%           H(s(n)) = [sqrt(x^2 + y^2 + z^2);
%                      acos(z/(sqrt(x^2 + y^2 + z^2))); 
%                      atan(y/x);
%                      (x*vx+y*vy+z*vz)/(sqrt(x^2 + y^2 + z^2))], and
%           v(n) is vector of measurement noise with covariance matrix
%           R(n)of size 6×6.
%
%   Step method syntax:
%
%   Y = step(H, U, UVAR, dt) updates tracker with measurment vector u =
%   [range azimuth velocity speed], variance vector uvar = [rangeVariance
%   azimuthVariance elevationVariance speedVariance], and time interval dt
%
%   radarTracking methods:
%
%   step                        - Tracking update
%
%   radarTracking properties:
%
%   maxNumOfTrackingObjects     - Maximum number of individual trackers
%   maxNumOfMeasurementsPerTracker - Maximum number of measurements that
%   can be associated with a single tracker
%   trackerActiveThreshold
%   trackerForgetThreshold
%   trackerAssociationThreshold - Maximum square distance between a measurment vector and associated tracker, in m^2
%   stateTrackingVectorType     - Supported types are '2D' or '3D' 

%   References
%   [1]	Moringa M., Shinoda H., and Kondoh H., “DOA Estimation of Coherent
%       Waves for 77 GHz Automotive Radar with Three receiving Antennas,”
%       Proceedings of the 6th European Radar Conference, 2009.
%   [2]	Kamen E. and Su. J.,  Introduction to Optimal Esitmation. Springer,
%       1999.

properties (Nontunable)
    %trackModuleConfig, received from upper layer
    trackerConfig
end

properties (Nontunable)
    sceneryNumBoundaryBoxes = 0;
    sceneryNumStaticBoxes = 0;
    
    %allocation, minimum total SNR of points in a set, ratio (default: -1, none)
    allocationMinTotalSNR = -1    
    %allocation, minimum total SNR of points in a set, ratio (default: -1, none)
    allocationMinTotalSNRObscured = -1    
    %allocation, minimal radial velocity of the target, m/s (default 0, none)
    allocationMinTargetVelocity = 0        
    %allocationThresholdPoints, minimum number of points (default: 3)
    allocationThresholdPoints = 5    
    %allocationThresholdDistance, in meters (default: 4)
    allocationThresholdDistance = 1
    %allocationThresholdVelocity, in m/c
    allocationThresholdVelocity = 2
    %unrollingRangeRateAlpha
    unrollingRangeRateAlpha = 0.5
    %unrollingRangeRateConf, in %
    unrollingRangeRateConf = 0.1
    %Axes
    ax
end

properties (Access = private, Nontunable)
    pHeartBeatCount
    pComputeFormat
    pStateVectorLength
    pMaxRadialVelocity
    pRadialVelocityResolution
    pStateVectorType
    pMeasurementVectorLength
    pMaxNumTracks
    pMaxNumPoints
    pMaxAcceleration
    pInitialRadialVelocity
    pProcessNoiseVariance
    pIgnoreTargetsBehindMe
    pMeasurementNoiseVariance
    pNumOfTrackingObjects
    pGatingMaxLimits
    pTrackingObjectList
    pAssociatedList
    pAllocatedList
    pRadialVelocity
    pTrackerIdList
    pPendingIndication
    
    pAllocationMinTotalSNR
    pAllocationMinTotalSNRObscured
    pAllocationMinTargetVelocity
    pAllocationThresholdPoints
    pAllocationThresholdDistance
    pAllocationThresholdVelocity
    
    pScenery
    
    pNumOfTotalTracks
    pColors
    pTracks
    pVerbose
    pDeltaT
    
    pUnitConfig
    pBestTrack
   
    hFigure
    hAxis2d
    hAxis3d    
    hMpoints
end


properties (Constant = true)
    STATE_FREE = 0;
    STATE_INIT = 1;
    STATE_DETECTION = 2;
    STATE_ACTIVE = 3;
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
    function obj = trackModule(varargin)
        setProperties(obj, nargin, varargin{:});

        % Create cell array to hold the track objects
        obj.pMaxNumTracks = obj.trackerConfig.maxNumTracks;
        obj.pMaxNumPoints = obj.trackerConfig.maxNumPoints;
        obj.pStateVectorType = obj.trackerConfig.stateVectorType;
        obj.pMaxRadialVelocity = obj.trackerConfig.maxRadialVelocity;
        obj.pRadialVelocityResolution = obj.trackerConfig.radialVelocityResolution;
        obj.pMaxAcceleration = obj.trackerConfig.maxAcceleration;
        obj.pInitialRadialVelocity = obj.trackerConfig.initialRadialVelocity;   

        obj.pDeltaT = obj.trackerConfig.deltaT;
        obj.hAxis2d = obj.trackerConfig.ax;
        
        switch(obj.trackerConfig.verbose)
            case 0
                % NONE
                obj.pVerbose = 0;
            case 1
                % ERRORS LEVEL
                obj.pVerbose = obj.VERBOSE_ERROR_INFO;
            case 2
                % WARNINGS LEVEL
                obj.pVerbose = sum([obj.VERBOSE_ERROR_INFO, obj.VERBOSE_WARNING_INFO]);
                obj.pVerbose = sum([obj.pVerbose,obj.VERBOSE_GATEXY_INFO, obj.VERBOSE_GATEG1_INFO]);
            case 3
                % DEBUG LEVEL
%                obj.pVerbose = sum([obj.VERBOSE_ERROR_INFO, obj.VERBOSE_WARNING_INFO, obj.VERBOSE_DEBUG_INFO, obj.VERBOSE_UNROLL_INFO, obj.VERBOSE_STATE_INFO]);
%                obj.pVerbose = sum([obj.pVerbose,obj.VERBOSE_GATEXY_INFO, obj.VERBOSE_GATEG1_INFO, obj.VERBOSE_ASSOSIATION_INFO]);
                obj.pVerbose = sum([obj.VERBOSE_STATE_INFO, obj.VERBOSE_GATEG1_INFO]);
%                 obj.pVerbose = obj.VERBOSE_STATE_INFO;
            case 4
                % ELEVATED LEVEL
                obj.pVerbose = sum([obj.VERBOSE_ERROR_INFO, obj.VERBOSE_WARNING_INFO, obj.VERBOSE_DEBUG_INFO, obj.VERBOSE_MATRIX_INFO]);
            case 5
                % MAXIMUM LEVEL
                obj.pVerbose = sum([obj.VERBOSE_ERROR_INFO, obj.VERBOSE_WARNING_INFO, obj.VERBOSE_DEBUG_INFO, obj.VERBOSE_MATRIX_INFO, obj.VERBOSE_UNROLL_INFO, obj.VERBOSE_STATE_INFO, obj.VERBOSE_ASSOSIATION_INFO]);
            case 6
                % VISUALIZATION
                obj.pVerbose = sum([obj.VERBOSE_ERROR_INFO, obj.VERBOSE_WARNING_INFO, obj.VERBOSE_DEBUG_INFO, obj.VERBOSE_MATRIX_INFO, obj.VERBOSE_UNROLL_INFO, obj.VERBOSE_STATE_INFO, obj.VERBOSE_ASSOSIATION_INFO]);
                obj.pVerbose = sum([obj.pVerbose,obj.VERBOSE_GATEXY_INFO, obj.VERBOSE_GATEG1_INFO]);  
        end  
                
        % Configure and Create all Track Units
        obj.pUnitConfig.maxNumPoints = obj.pMaxNumPoints;
        obj.pUnitConfig.stateVectorType = obj.pStateVectorType;
        obj.pUnitConfig.maxRadialVelocity = obj.pMaxRadialVelocity;
        obj.pUnitConfig.radialVelocityResolution = obj.pRadialVelocityResolution;
        obj.pUnitConfig.maxAcceleration = obj.pMaxAcceleration;
        obj.pUnitConfig.initialRadialVelocity = obj.pInitialRadialVelocity;
        obj.pUnitConfig.verbose = obj.pVerbose;
        obj.pUnitConfig.deltaT = obj.pDeltaT;
        
        if(~isempty(obj.trackerConfig.advParams.allocation))
            obj.pAllocationMinTotalSNR = obj.trackerConfig.advParams.allocation(1);
            obj.pAllocationMinTotalSNRObscured = obj.trackerConfig.advParams.allocation(2);
            obj.pAllocationMinTargetVelocity = obj.trackerConfig.advParams.allocation(3);
            obj.pAllocationThresholdPoints = obj.trackerConfig.advParams.allocation(4);
            obj.pAllocationThresholdDistance = obj.trackerConfig.advParams.allocation(5);
            obj.pAllocationThresholdVelocity = obj.trackerConfig.advParams.allocation(6);
        else
            obj.pAllocationMinTotalSNR = obj.allocationMinTotalSNR;
            obj.pAllocationMinTotalSNRObscured = obj.allocationMinTotalSNRObscured;
            obj.pAllocationMinTargetVelocity = obj.allocationMinTargetVelocity;
            obj.pAllocationThresholdPoints = obj.allocationThresholdPoints;
            obj.pAllocationThresholdDistance = obj.allocationThresholdDistance;
            obj.pAllocationThresholdVelocity = obj.allocationThresholdVelocity;            
        end
        
        % scenery is always sent from module to unit
        if(~isempty(obj.trackerConfig.advParams.scenery))
            ind = 1;
            obj.pScenery.numBoundaryBoxes = obj.trackerConfig.advParams.scenery(ind);
            ind = ind + 1;
            if(obj.pScenery.numBoundaryBoxes ~=0)
                obj.pScenery.boundaryBoxes = reshape(obj.trackerConfig.advParams.scenery(ind:ind+obj.pScenery.numBoundaryBoxes*4-1),obj.pScenery.numBoundaryBoxes,4);
            end
            ind = ind + obj.pScenery.numBoundaryBoxes*4;
            obj.pScenery.numStaticBoxes = obj.trackerConfig.advParams.scenery(ind);
            ind = ind + 1;
            if(obj.pScenery.numStaticBoxes ~=0)
                obj.pScenery.staticBoxes = reshape(obj.trackerConfig.advParams.scenery(ind:ind+obj.pScenery.numBoundaryBoxes*4-1),obj.pScenery.numBoundaryBoxes,4);
            end           
            obj.pUnitConfig.advParams.scenery = obj.pScenery;
        else
            obj.pScenery.numBoundaryBoxes = obj.sceneryNumBoundaryBoxes;
            obj.pScenery.numStaticBoxes = obj.sceneryNumStaticBoxes;
            obj.pUnitConfig.advParams.scenery = obj.pScenery;            
        end
        
        if(~isempty(obj.trackerConfig.advParams.gating))
            obj.pUnitConfig.advParams.gating = obj.trackerConfig.advParams.gating;
        else
            obj.pUnitConfig.advParams.gating = []; % Use defaults
        end
        if(~isempty(obj.trackerConfig.advParams.thresholds))
            obj.pUnitConfig.advParams.thresholds = obj.trackerConfig.advParams.thresholds;
        else
            obj.pUnitConfig.advParams.thresholds = []; % Use defaults
        end        
        if(~isempty(obj.trackerConfig.advParams.variations))
            obj.pUnitConfig.advParams.variations = obj.trackerConfig.advParams.variations;
        else
            obj.pUnitConfig.advParams.variations = []; % Use defaults
        end
        
        % Create all tracking units
        obj.pTracks = cell(1,obj.pMaxNumTracks);
        for tid = 1:obj.pMaxNumTracks
            obj.pUnitConfig.tid = tid;
            obj.pTracks{tid} = trackUnit('trackUnitConfig', obj.pUnitConfig);
        end
        
        if(bitand(obj.pVerbose, obj.VERBOSE_GATERA_INFO))
            obj.hFigure = figure('Name','Tracking Gating and Association','NumberTitle','off');
            axis([0 80, -0.1 1, -50 0]);
            grid on; hidden off; hold on;
            obj.hAxis3d = obj.hFigure.CurrentAxes;  
        end
        
        obj.hMpoints = 0;
    end
    
    function handle = getTrackerhandle(obj, id)
        handle = obj.pTracks{id};
    end
    
    function setTrackerMark(obj, id, value)
        trackUnitSetMarker(obj.pTracks{id}, value);
    end    
    function marker = getTrackerMark(obj, id)
        marker = trackUnitGetMarker(obj.pTracks{id});
    end    
end
methods (Access = protected)
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@matlab.System(obj);
    end
    
    function loadObjectImpl(obj,s,~)
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end
    
    function resetImpl(obj)
        obj.pNumOfTrackingObjects = 0;
        obj.pNumOfTotalTracks = 0;
        obj.pTrackingObjectList = [];
        obj.pTrackerIdList = 1:obj.pMaxNumTracks;
    end
    
    function setupImpl(obj)        
        obj.pHeartBeatCount = 0;
        obj.pColors = 'bgrcm';
    end

    function [track, numOfReportedTracks, mIndex] = stepImpl(obj, u, uvar, mNum)

        if(mNum == 0)
            u = [];
            uvar = [];
        end
        
        obj.pHeartBeatCount = obj.pHeartBeatCount + 1;

        if(~isempty(obj.hAxis2d))
            cla(obj.hAxis2d);
            hold(obj.hAxis2d, 'on');
        end
        if(~isempty(obj.hAxis3d))
            cla(obj.hAxis3d);
        end
        
        % For each step, we will re-associate measurements with trackers. 
        obj.pAssociatedList = zeros(obj.pMaxNumTracks,obj.pMaxNumPoints);
        obj.pAllocatedList = zeros(obj.pMaxNumPoints, 1);
        obj.pRadialVelocity = zeros(obj.pMaxNumPoints, 1);
        
        mIndex = ones(1,size(u,2))*255;

        % Ignore points outside boundary boxes
        if(mNum)
            posAll = [u(1,:).*sin(u(2,:)); u(1,:).*cos(u(2,:))];
            inBoxInd = false(1,size(posAll,2));
            for n = 1:obj.pScenery.numBoundaryBoxes
                box = obj.pScenery.boundaryBoxes(n,:);
                inBoxInd = or(inBoxInd, (posAll(1,:) > box(1)) & (posAll(1,:) < box(2)) & (posAll(2,:) > box(3)) & (posAll(2,:) < box(4)));            
            end
            u = u(:,inBoxInd);        
        end
        
        % Create measurements update pending indications
        obj.pPendingIndication = ones(1, size(u,2));

        predictTrackers(obj);        
        associateTrackers(obj, u);
        allocateTrackers(obj, u);  
        updateTrackers(obj, u, uvar);
        [track, numOfReportedTracks] = reportTrackers(obj);
        if(mNum)
            mIndex(inBoxInd) = obj.pBestTrack;
        end    
%        text(obj.hAxis2d, 0,76, ['Total #', num2str(obj.pNumOfTotalTracks), ' targets']);
%        text(obj.hAxis2d, 0,74, ['Tracking #', num2str(numOfReportedTracks), ' targets']);
    end
    
    function predictTrackers(obj)
        for n = 1:length(obj.pTrackingObjectList)
            tid = obj.pTrackingObjectList(n);            
            trackUnitPredict(obj.pTracks{tid});
        end
    end
   
    function associateTrackers(obj, u)
        nMeas = size(u,2);
        if(nMeas)
            bestTrack = zeros(nMeas,1);
            % Visualize measurment points
            if(bitand(obj.pVerbose, obj.VERBOSE_ASSOSIATION_INFO))
                if(~isempty(obj.hAxis2d))
                    plot3(obj.hAxis2d, u(1,:).*sin(u(2,:)), u(1,:).*cos(u(2,:)), u(3,:), '.k');
                end
                if(~isempty(obj.hAxis3d))
                    plot3(obj.hAxis3d, u(1,:), u(2,:), u(3,:), '.k');
                end
            end
        end
        
        nTrack = length(obj.pTrackingObjectList);
        if((nMeas ~= 0) && (nTrack ~= 0))
            % We have both: active trackers AND pending measurments
            bestScore = inf(nMeas,1);
            
            for n=1:nTrack
                
                tid = obj.pTrackingObjectList(n);
                [bestScore, bestTrack, rvOut, centroid, ec, G] = trackUnitAssociate(obj.pTracks{tid}, u, bestScore, bestTrack);
                ind = (bestTrack==tid);
                obj.pRadialVelocity(ind) = rvOut(ind);
                
                % Visualize Track centroid
                color = obj.pColors(mod(tid,length(obj.pColors))+1);
                if(bitand(obj.pVerbose, obj.VERBOSE_ASSOSIATION_INFO))
                    if(~isempty(obj.hAxis2d))
                        plot3(obj.hAxis2d, centroid(1)*sin(centroid(2)),centroid(1)*cos(centroid(2)),centroid(3), '*', 'color', color);                 
                    end
                    if(~isempty(obj.hAxis3d))
                        plot3(obj.hAxis3d, centroid(1)*sin(centroid(2)),centroid(1)*cos(centroid(2)),centroid(3), '*', 'color', color);                 
                        plot3(obj.hAxis3d, centroid(1),centroid(2),centroid(1), '*', 'color', color);                 
                    end
                end                
                if(bitand(obj.pVerbose, obj.VERBOSE_GATEXY_INFO) || bitand(obj.pVerbose, obj.VERBOSE_GATERA_INFO) || bitand(obj.pVerbose, obj.VERBOSE_GATEG1_INFO))
                    v = computeVolume(ec);
                    gUnitary = (1/v)^(2/3);
                    [x1, y1, z1, v1] = gatePlot3(obj.hAxis2d, gUnitary, centroid, ec);
                    [xG, yG, zG, vG] = gatePlot3(obj.hAxis2d, G, centroid, ec);
                end
                
                if(bitand(obj.pVerbose, obj.VERBOSE_GATEG1_INFO))
                    if(~isempty(obj.hAxis2d))
                        h1 = mesh(obj.hAxis2d, x1.*sin(y1),x1.*cos(y1),z1);
                        h1.FaceColor = 'none';
                        h1.EdgeColor = [0.5 0.5 0.5];
                    end
                end
                
                if(bitand(obj.pVerbose, obj.VERBOSE_GATEXY_INFO))                    
                    if(~isempty(obj.hAxis2d))
                        h2 = mesh(obj.hAxis2d, xG.*sin(yG),xG.*cos(yG),zG);
                        h2.FaceColor = 'none';
                        h2.EdgeColor = color;
                        view(obj.hAxis2d, [0,90]); 
                    end
                end
                
                if(bitand(obj.pVerbose, obj.VERBOSE_GATERA_INFO))
                    if(~isempty(obj.hAxis3d))
                        h1 = mesh(obj.hAxis3d, x1, y1, z1);
                        h1.FaceColor = 'none';

                        h2 = mesh(obj.hAxis3d, xG, yG, zG);
                        h2.FaceColor = 'none';
                        view(obj.hAxis3d, [0,90]);
                    end
                end                
            end
%{            
            % Associate one measurement at a time to the closest track
            % Multiple measurements can be associated to one track
            % A measurement may only be associated to one and only one track. 
            % A measurement can be left non-assiciated
            [mScore, closestTracker2] = min(score,[],1);
            for mid=1:length(closestTracker2)
                val = mScore(mid);
                tid = closestTracker2(mid);
%                if(val < obj.pTrackerAssociationThreshold)
                if(val < 100)
                    id = obj.pTrackingObjectList(tid);
                    ind = find(obj.pAssociatedList(id,:)==0,1);
                    obj.pAssociatedList(id,ind) = mid;
                    obj.pPendingIndication(mid) = 0;
                end
            end
 %}
            for id=1:nTrack
                tid = obj.pTrackingObjectList(id);
                mid = find(bestTrack==tid);
                obj.pAssociatedList(tid,1:length(mid))=mid;
                obj.pPendingIndication(mid) = 0;
            end
                
            % Visualize association process
            for id=1:nTrack
                tid = obj.pTrackingObjectList(id);
                aid = obj.pAssociatedList(tid,:)~=0;
                mid = obj.pAssociatedList(tid,aid);                
                color = obj.pColors(mod(tid,length(obj.pColors))+1);
                if(bitand(obj.pVerbose, obj.VERBOSE_ASSOSIATION_INFO))
                    if(~isempty(obj.hAxis2d))                        
                        plot3(obj.hAxis2d, u(1,mid).*sin(u(2,mid)), u(1,mid).*cos(u(2,mid)),u(3,mid), 'o','color', color);
                    end
                    if(~isempty(obj.hAxis3d))
                        plot3(obj.hAxis3d, u(1,mid),u(2,mid),u(3,mid), 'o','color', color);
                    end
                end
            end
            
            % How many points left unassociated
            numLeft = nnz(obj.pPendingIndication);
            if(numLeft)
                if(bitand(obj.pVerbose, obj.VERBOSE_WARNING_INFO))                
                    disp(['Warning: ', num2str(numLeft), ' measurements left un-associated']);
                end
            end
        end
        
        if(nMeas)
            obj.pBestTrack = bestTrack;
        else
            obj.pBestTrack = [];
        end
    end
    
    function allocateTrackers(obj, u)
       totalSNR = 0;
       for n = 1:length(obj.pPendingIndication)            
            if(obj.pPendingIndication(n) == 1)
           
                obj.pPendingIndication(n) = 0;
                obj.pAllocatedList = zeros(obj.pMaxNumPoints, 1);
                obj.pAllocatedList(1) = n;
                numPoints = 1;
                if(size(u,1) > 3)
                    totalSNR = u(4,n);
                end
                um = u(1:3,n);

                for k = n+1:length(obj.pPendingIndication)
                    if(obj.pPendingIndication(k) == 1)
                        rv = unrollRadialVelocity(obj.pMaxRadialVelocity, um(3), u(3,k));
                        if(abs(rv - um(3)) < obj.pAllocationThresholdVelocity)                            
                            dist = u(1,k)^2 + um(1)^2 - 2*u(1,k)*um(1)*cos(um(2)-u(2,k));
                            if(sqrt(dist) < obj.pAllocationThresholdDistance)
                                ind = find(obj.pAllocatedList==0,1); %find first available slot
                                obj.pAllocatedList(ind) = k;
%                                obj.pPendingIndication(k) = 0;
                                
                                ind = obj.pAllocatedList~=0;
                                mid = obj.pAllocatedList(ind);
                                um = mean(u(:,mid),2);
%                               um = (max(u(:,m),[],2)+min(u(:,m),[],2))/2;
                                numPoints = numPoints + 1;
                                if(size(u,1) > 3)
                                    totalSNR = totalSNR + u(4,k);
                                end                                
                            end
                        end
                    end
                end
                
                if( (numPoints > obj.pAllocationThresholdPoints) && ...
                    (abs(um(3)) > obj.pAllocationMinTargetVelocity) )

                    %Ghost behind the target protection
                    passSNRThreshold = false;
                    isObscured = false;
                    for id = 1:length(obj.pTrackingObjectList)
                        tid = obj.pTrackingObjectList(id);
                        if(trackUnitIsBehind(obj.pTracks{tid}, um))
                            isObscured = true;
                            break;
                        end
                    end
                    
                    if isObscured
                        if(totalSNR > obj.pAllocationMinTotalSNRObscured)
                           passSNRThreshold = true;
                        end
                    else
                        if(totalSNR > obj.pAllocationMinTotalSNR)
                           passSNRThreshold = true;
                        end
                    end
                    % Allocate only if SNR threshold passed
                    if passSNRThreshold    
                        tid = allocateTrackerID(obj); 
                        trackUnitInitialize(obj.pTracks{tid}, obj.pHeartBeatCount, um);

                        obj.pAssociatedList(tid,:) = obj.pAllocatedList;
                        obj.pPendingIndication(mid) = 0;
                    end
                end                    
            end
       end
    end

    function updateTrackers(obj, u, uvar)
        trackingList = obj.pTrackingObjectList;
                
        for n=1:length(trackingList)
            tid = trackingList(n);
            
            if(obj.pAssociatedList(tid,1) ~= 0)
                % Tracker is associated with at least one measurement              
                ind = (obj.pAssociatedList(tid,:) ~= 0);
                mid = obj.pAssociatedList(tid,ind);
                u(3,mid) = obj.pRadialVelocity(mid);
            else
                mid = [];
            end

            state = trackUnitUpdate(obj.pTracks{tid}, u, uvar, mid);  
            
            if(state == obj.STATE_FREE)
                obj.pNumOfTotalTracks = obj.pNumOfTotalTracks +1;
                if(bitand(obj.pVerbose, obj.VERBOSE_DEBUG_INFO))
                    disp(['Tick ',num2str(obj.pHeartBeatCount),': tracker #', num2str(tid),' freed', ]);
                end
                freeTrackerID(obj, tid);
            end
        end
    end
    
    function [track, numOfReportedTracks] = reportTrackers(obj)
        trackStruct = struct('S',zeros(1,obj.pStateVectorLength), 'trackerId', 0, 'measureId', 0, 'ec', zeros(obj.pMeasurementVectorLength, obj.pMeasurementVectorLength));
        track = repmat(trackStruct, 1, obj.pMaxNumTracks);
        numOfReportedTracks = 0;
        for n=1:length(obj.pTrackingObjectList)
            tid = obj.pTrackingObjectList(n);
            numOfReportedTracks = numOfReportedTracks + 1;
            [S, state, ec, g] = trackUnitGetInfo(obj.pTracks{tid});
            track(numOfReportedTracks).state = state;
            track(numOfReportedTracks).S = S;
            track(numOfReportedTracks).trackerId = tid;
            track(numOfReportedTracks).ec = ec;            
            track(numOfReportedTracks).g = g;            
            ind = (obj.pAssociatedList(tid,:) ~= 0);
            mid = obj.pAssociatedList(tid,ind);
            track(numOfReportedTracks).measureId = mid;
        end
    end

    function id = allocateTrackerID(obj)
    %allocateTrackerID Allocate new tracker ID from the free pool
        % Get new tracker ID
        [id, obj.pTrackerIdList] = getFromList(obj.pTrackerIdList);
        obj.pNumOfTrackingObjects = obj.pNumOfTrackingObjects + 1;
        % Adds tracker to the list
        obj.pTrackingObjectList = addToList(obj.pTrackingObjectList, id);
    end

    function freeTrackerID(obj, id)
    %freeTrackerID Places tracker ID back to free pool
        % Remove tracker from the list
        obj.pTrackingObjectList = removeFromList(obj.pTrackingObjectList, id);
        % Returns tracker id back to the pool
        obj.pTrackerIdList = addToList(obj.pTrackerIdList, id);
        obj.pNumOfTrackingObjects = obj.pNumOfTrackingObjects - 1;
    end
end
end

% List functions
function list = addToList (list, id)
    list = [list id];
end   
function list = removeFromList (list, id)
    ind = (list==id);
    list(ind) = [];    
end
function [id, list] = getFromList (list)
    if(isempty(list))
        id = 0;
    else
        id = list(1);
        list(1) = [];
    end
end
function v = computeVolume(A)
    a = 1/sqrt(A(1,1));
    b = 1/sqrt(A(2,2));
    c = 1/sqrt(A(3,3));
    v = 4*pi/3*a*b*c;
end
