close all;
clear all;
colors='bgrcm';

diary('tmLog.txt');

simMinutesRun = 10; %60 One hour of logging
simPlayData = 0;
capturedData = 0;
if(capturedData == 1)
    
    load('aoaInfoDs2cfar1.mat');
    n =1;
    while(isempty(aoaInfo{n}))
        n = n+1;
    end
    capturedFrameNum = n;
    capturedEndFrame = 500;
end

tDataStruct = struct('aInd', [], 'vType',[], 'vSpeed', []);
if (simPlayData==1)   
    load('testTracking-Split000.mat', 'tData');
%    load('testTracking-RangeSeparation006.mat', 'tData');    
%    load('testTracking-VelSeparation000.mat', 'tData');    
else
    tData = repmat(tDataStruct,simMinutesRun, 1);
end

vHistStruct = struct('veh', [], 'num', []);
vHist = repmat(vHistStruct,4, 1);
tHistStruct = struct('tid', 0, 'vid', 0, 'allocationTime', 0, 'tick', 0, 'sHat', '[]');
tHist_current = repmat(tHistStruct,1, 100);

pCloudStruct = struct('num', [], 'u', [], 'uvar', [], 'vid', []); 
cHist = repmat(pCloudStruct, 1, 1500);

sceneName = 'tmTracking';
videoFrame = struct('cdata',[],'colormap', []);

load('randomState2.mat', 'randomState');
rng(randomState);
%randomState = rng;
%save('randomState.mat');

laneStatStruct = struct('total',[], 'moving',[], 'pending', [], 'pcloud',[]);

%sensor parameters
sensor.rangeMax = 70;
sensor.rangeMin = 15;
sensor.azimuthFoV = 30*pi/180; %30 degree FOV in horizontal direction
sensor.azimuthTilt = asin(10/70);
sensor.azimuthMin = (pi/2 - (sensor.azimuthTilt + sensor.azimuthFoV/2));
sensor.azimuthMax = sensor.azimuthMin+sensor.azimuthFoV;
sensor.framePeriod = 40;
sensor.maxURadialVelocity = 20;

%scene parameters
scene.numberOfLanes = 4;
scene.laneWidth = [3 3 3 3];
scene.laneXPos = [-3 0 3 6];
scene.stopLineDist = [21 20 20 20];
scene.stopExitDist = [15 15 15 15];
scene.trafficTimeouts = [20 2 5]; % [20 2 10] Green time, Yellow time, Red Time, (s)
%scene.lambda = [1/10, 2/10, 2/10, 3/10]; % arrival rate per second, (1/10 => 1 vehicle in 10 seconds)
scene.lambda = [1/10, 2/10, 3/10, 4/10]; % arrival rate per second, (1 vehicle in 10 seconds)
%scene.vehicleProb = [0.2 0.2 0.1 0.5; 0.3 0.6 0.1 0; 0.5 0.4 0.1 0; 0.9 0 0.1 0]; %Cars/Trucks/Motorcycles/Bicycles
scene.vehicleProb = [0.5 0 0 0.5; 1 0 0 0; 1 0 0 0; 1 0 0 0]; %Cars/Trucks/Motorcycles/Bicycles

dt = sensor.framePeriod/1000;
hPlotCloudHandle = 0;

scrsz = get(groot,'ScreenSize');
hFigure = figure('Name','Traffic Monitoring','NumberTitle','off','Position',[0 0 1600 900]);

hStatGlobal = zeros(1,3);
hRb = uicontrol(hFigure,'style','radio','string','Pause',...
                 'Position',[10 50 50 20],...
                 'Value',1,'UserData',1);

angle=linspace(sensor.azimuthMin, sensor.azimuthMax, 128);
for n=1:4
    subplot(1,4,n);
    plot([0 sensor.rangeMax*cos(angle) 0],[0 sensor.rangeMax*sin(angle) 0], '-k'); hold on
    plot(sensor.rangeMin*cos(angle),sensor.rangeMin*sin(angle), '-k');    
    xlabel('X coordinate, m');
    ylabel('Y coordinate, m');
    %axis([-rangeMax*1.2 rangeMax*1.2 0 rangeMax*1.2]);
    axis([-10 30 0 80]);
    for nLane = 1:scene.numberOfLanes
        if(nLane == 1)
            line('Xdata',[scene.laneXPos(nLane) scene.laneXPos(nLane)],'YData',[scene.stopLineDist(nLane) 80], 'Color','k', 'LineStyle', '-', 'LineWidth', 1);
        else
            line('Xdata',[scene.laneXPos(nLane) scene.laneXPos(nLane)],'YData',[scene.stopLineDist(nLane) 80], 'Color','k', 'LineStyle', '--');
        end
        line('Xdata',[scene.laneXPos(nLane) scene.laneXPos(nLane)+scene.laneWidth(nLane)],'YData', [scene.stopLineDist(nLane) scene.stopLineDist(nLane)],'Color','r');            
    end
    line('Xdata',[scene.laneXPos(nLane)+scene.laneWidth(nLane) scene.laneXPos(nLane)+scene.laneWidth(nLane)],'YData',[scene.stopLineDist(nLane) 80], 'Color','k', 'LineStyle', '-', 'LineWidth', 1);
    grid on;
    grid minor;
    ax(n) = gca;
end

% Initialize traffic light object

tmLight = trafficLight('timeouts', scene.trafficTimeouts);
trafficLightColors = 'gyr';
trafficLights = {'Green', 'Yellow', 'Red'};
%trafficLightPos = [-6 18.5; -6 20; -6 21.5];
trafficLightPos = [-6 18; -6 20; -6 22];
hLight(1) = plot(ax(1), trafficLightPos(1,1),trafficLightPos(1,2),'o','MarkerSize',10,'MarkerFaceColor', 'g', 'MarkerEdgeColor','k');
hLight(2) = plot(ax(1), trafficLightPos(2,1),trafficLightPos(2,2),'o','MarkerSize',10,'MarkerFaceColor', 'w', 'MarkerEdgeColor','k');
hLight(3) = plot(ax(1), trafficLightPos(3,1),trafficLightPos(3,2),'o','MarkerSize',10,'MarkerFaceColor', 'w', 'MarkerEdgeColor','k');
activeLight = 1;
vehicleType{scene.numberOfLanes} = 0;
hStatLane = zeros(2, scene.numberOfLanes);
vehicleCount = int64(zeros(scene.numberOfLanes,1));
vehicleCountArrived = int64(zeros(scene.numberOfLanes,1));
objectCount = int64(zeros(scene.numberOfLanes,1));


% Initialize pointCloud object
%load('cloudStatMax.mat', 'cloudStat'); % Average number of detections as a function of the range, collected experimentally
load('cloudStat-cfar2-dw-02.mat', 'cloudStat'); % Average number of detections as a function of the range, collected experimentally
cloud.radarPos = [0,0];
cloud.uRadVelocityMax = sensor.maxURadialVelocity; % Unambigious Radial Velocity (+/- 10m/s), for testing purposes
cloud.dopplerFftBinRes = 0.3298;
cloud.stat.ranges = cloudStat.ranges;
cloud.stat.numPoints = cloudStat.numPoints;
cloud.stat.xStd = cloudStat.xStd;
cloud.stat.ycdf = cloudStat.ycdf;
cloud.stat.yval = cloudStat.yval;
tmCloud = pointCloud('cloud', cloud);

% Initialize lane objects
tmLane = cell(scene.numberOfLanes,1);
for nLane = 1:scene.numberOfLanes
    lane.num = nLane;
    lane.width = scene.laneWidth(nLane);
    lane.xpos = scene.laneXPos(nLane);
    lane.enterline = sqrt(sensor.rangeMax^2-lane.xpos^2);
    lane.stopline = scene.stopLineDist(nLane);
    lane.exitline = scene.stopExitDist(nLane);
    % create traffic lane objects
    tmLane{nLane} = trafficLane('lane', lane, 'pCloudHandle', tmCloud, 'ax', ax);
end

% Create tracking module
maxNumTracks = 100;
maxNumPoints = 1000;

trackerConfig.stateVectorType = 1; % 0=>2D, 1=>2DA
trackerConfig.verbose = 2;
trackerConfig.maxNumPoints = maxNumPoints;
trackerConfig.maxNumTracks = maxNumTracks;

trackerConfig.maxURadialVelocity = sensor.maxURadialVelocity;
trackerConfig.maxAcceleration = 20;
trackerConfig.deltaT = sensor.framePeriod/1000;
trackerConfig.initialRadialVelocity = -20;
trackerConfig.ax = ax(4);

trackerConfig.advParams.allocation = [-1, 0, 3, 4, 2]; % any SNR, 1m/s doppler, a set of 3 Points with 4m in distance and 2m/s velocity separation
trackerConfig.advParams.gating = [16,12,8,0]; % 16 units of Volume, Limits: 12m in length, 8m in width and no limit it velocity
trackerConfig.advParams.thresholds = [3,3, 5,5,1]; % det2active, det2forget, act2free, stat2free, exit2free
trackerConfig.advParams.variations = [1/sqrt(12), 1/sqrt(12), 2]; % Height, Width, and Doppler std
trackerConfig.advParams.scenery = [-100,100,20,100]; % No walls, 10m lower exit zone, no upper zone. 

hRadarTracking = trackModule('trackerConfig', trackerConfig);

hPlotTrackHandle = zeros(maxNumTracks, 1);
trackedPos{maxNumTracks} = [];
objectTrackingState = zeros(maxNumTracks, 1);
activeTracks = zeros(1, trackerConfig.maxNumTracks);

%hCloudSubplot = subplot(1,2,2);
title(ax(1), 'Simulated Scene');
title(ax(2), 'Simulated Point Cloud');
title(ax(3), 'Tracking History');
title(ax(4), 'Tracking Gating and Association');


for minuteTick=1:simMinutesRun
    disp(minuteTick);
    F = repmat(videoFrame,1500,1);
    
    if simPlayData==1
        aInd = tData(minuteTick).aInd;
        vType = tData(minuteTick).vType;
        vSpeed = tData(minuteTick).vSpeed;
    else
        % generate 60 seconds of vehicle arrivials data
        aInd = genArrivals(60, dt, scene.lambda);
        if(minuteTick == simMinutesRun)
            %If that is the last minute, then zero down last 100 ticks
            %arrivals to allow all objects to propogate
            aInd(end-100:end,:) = 0;
        end
        vType = genVehicleTypes(sum(aInd), scene.vehicleProb);
        vSpeed = cell(scene.numberOfLanes,1);
        for n=1:scene.numberOfLanes
            vSpeed{n}=zeros(1,nnz(aInd(:,n)));
        end
        tData(minuteTick).aInd = aInd;
        tData(minuteTick).vType = vType;
        tData(minuteTick).vSpeed = vSpeed;
    end

    nArrivals = ones(scene.numberOfLanes,1);
    stat = repmat(laneStatStruct,scene.numberOfLanes, 1);
  
    for tick = 1:size(aInd,1)
        gtick = (minuteTick-1)*size(aInd,1)+tick;

        if(gtick == 670)
            disp(gtick);
        end
%        disp(gtick);

        pCloudTm.id = [];       
        pCloudTm.xy = [];
        pCloudTm.u = [];
        pCloudTm.uvar = [];

        if(capturedData)
%            if(aoaInfo{capturedFrameNum}.numDynTargets)
%                ind = aoaInfo{capturedFrameNum}.confDopplerVal ~=0;
                pCloudTm.u = [angleEst{capturedFrameNum}.range; angleEst{capturedFrameNum}.angles*pi/180; angleEst{capturedFrameNum}.doppler];
                pCloudTm.xy = [pCloudTm.u(1,:).*sin(pCloudTm.u(2,:)); pCloudTm.u(1,:).*cos(pCloudTm.u(2,:))];
                pCloudTm.uvar = [1,(pi/180)^2,0.25]'.*ones(size(pCloudTm.u));
                %pCloudTm.uvar = (pCloudTm.u - mean(pCloudTm.u, 2)).^2;
 %           end
            capturedFrameNum = capturedFrameNum +1;
            if(capturedFrameNum > capturedEndFrame)
                break;
            end
        else
            for nLane = 1:scene.numberOfLanes
                arrival.tick = gtick;
                if(aInd(tick,nLane))
                    switch vType{nLane}(nArrivals(nLane))
                        case 1
                            if(vSpeed{nLane}(nArrivals(nLane)) == 0)
                                arrival.speed = ceil(randn(1,1)*10+25);
                            else
                                arrival.speed = vSpeed{nLane}(nArrivals(nLane));
                            end
                            arrival.type = 'car';
                        case 2
                            if(vSpeed{nLane}(nArrivals(nLane)) == 0)
                                arrival.speed = ceil(randn(1,1)*10+20);
                            else
                                arrival.speed = vSpeed{nLane}(nArrivals(nLane));
                            end
                            arrival.type = 'truck';
                        case 3
                            if(vSpeed{nLane}(nArrivals(nLane)) == 0)
                                arrival.speed = ceil(randn(1,1)*10+25);
                            else
                                arrival.speed = vSpeed{nLane}(nArrivals(nLane));
                            end
                            arrival.type = 'motorcycle';
                        case 4
                            if(vSpeed{nLane}(nArrivals(nLane)) == 0)
                                arrival.speed = ceil(randn(1,1)*3+10);
                            else
                                arrival.speed = vSpeed{nLane}(nArrivals(nLane));
                            end
                            arrival.type = 'bicycle';
                    end
                    vehicleCountArrived(nLane) = vehicleCountArrived(nLane) + 1;
                    nArrivals(nLane) = nArrivals(nLane)+1;
                else
                    arrival.type = 'none';
                    arrival.speed = 0;
                end
                
                % Process traffic in each Lane
                stat(nLane) = step(tmLane{nLane}, arrival, trafficLights{activeLight}, dt);
                if(hStatLane(1,nLane) ~= 0)
                    delete(hStatLane(1,nLane));
                end
                if(hStatLane(2,nLane) ~= 0)
                    delete(hStatLane(2,nLane));
                end
                hStatLane(1,nLane) = text(ax(1), scene.laneXPos(nLane)+0.4,74, [num2str(stat(nLane).moving), '/', num2str(stat(nLane).total)]);
                hStatLane(2,nLane) = text(ax(1), scene.laneXPos(nLane)+0.2,76, num2str(vehicleCount(nLane),'%03d'));

                if ~isempty(stat(nLane).pcloud.u)
                    inRangeInd = (stat(nLane).pcloud.u(1,:) > sensor.rangeMin) & (stat(nLane).pcloud.u(1,:) < sensor.rangeMax);
                
                    pCloudTm.id = [pCloudTm.id, stat(nLane).pcloud.id(:,inRangeInd)];
                    pCloudTm.xy = [pCloudTm.xy, stat(nLane).pcloud.xy(:,inRangeInd)];
                    pCloudTm.u = [pCloudTm.u, stat(nLane).pcloud.u(:,inRangeInd)];
                    pCloudTm.uvar = [pCloudTm.uvar, stat(nLane).pcloud.uvar(:,inRangeInd)];
                end   
            end
        end
        
        if(hPlotCloudHandle ~=0)
            delete(hPlotCloudHandle);
            hPlotCloudHandle = 0;
        end
        
        if ~isempty(pCloudTm.xy)
            hPlotCloudHandle = plot(ax(3), pCloudTm.xy(1,:), pCloudTm.xy(2,:),'.k');
        end
        
        numOutputPoints = size(pCloudTm.u,2);
        pCloudTm.u = [pCloudTm.u; zeros(1,numOutputPoints)]; % Add SNRs

        cHist(tick).num = numOutputPoints;
        cHist(tick).u = pCloudTm.u;
        cHist(tick).uvar = pCloudTm.uvar;
        cHist(tick).vid = pCloudTm.id;
        
%        [track, num] = step(hRadarTracking, pointCloud, pCloudTm.id);
        
        [track, num] = step(hRadarTracking, pCloudTm.u, pCloudTm.uvar, numOutputPoints);        
        
        trackedIDs = zeros(maxNumTracks, 1);

        % Collect tracking hidtory
        if(num)
            for m=1:num
                tid = track(m).trackerId;
                if(activeTracks(tid) == 0)
                    activeTracks(tid) = 1; % mark as reported
                    tHist_current(tid).tid = tid;
 %                   tHist_current(tid).vid = track(m).vid;
                    tHist_current(tid).allocationTime = gtick;
                    tHist_current(tid).tick = 1;
                    tHist_current(tid).sHat = zeros(1000,6);
                    tHist_current(tid).sHat(1,:) = track(m).S;
                else
                    activeTracks(tid) = 1; % mark as reported
                    ind = tHist_current(tid).tick + 1;
                    tHist_current(tid).sHat(ind,:) = track(m).S;
                    tHist_current(tid).tick = ind;
                end
            end
            ind = (activeTracks == 2); % not reported => ready to safe
            activeTracks(ind) = 3;
            ind = (activeTracks == 1); % all reported => active
            activeTracks(ind) = 2;
        end        
        
        if(num)
            % Display tracking trajectory
            for m=1:num
                tid = track(m).trackerId;
                trackedIDs(tid) = 1;
                mid = track(m).measureId;
                if(objectTrackingState(tid) == 0)
                    if(track(m).S(2) > 30)
                        objectTrackingState(tid) = 1;
                    end
                end
                if(track(m).S(2) < 25)
                    if (objectTrackingState(tid) == 1)
                        % count the object                     
                        if((track(m).S(1) > scene.laneXPos(1)) && (track(m).S(1) < (scene.laneXPos(end) + scene.laneWidth(end))))
                            ind = find(scene.laneXPos < track(m).S(1),1, 'last');
                            objectCount(ind) = objectCount(ind) +1;
                        end
                        objectTrackingState(tid) = 2;
                    end
                end
                % update trajectory
                trackedPos{tid} = [trackedPos{tid}; track(m).S(1), track(m).S(2)];
                if(hPlotTrackHandle(tid)~=0)
                    delete(hPlotTrackHandle(tid));
                end
                hPlotTrackHandle(tid) = plot(ax(3), trackedPos{tid}(:,1), trackedPos{tid}(:,2), '.', 'color', colors(mod(tid,length(colors))+1));
            end
        end
        ind = (trackedIDs == 0);
        objectTrackingState(ind) = 0;
        
        % Delete lost track and its trajectory
        ind = (hPlotTrackHandle~=0) & (trackedIDs == 0);        
        if(hPlotTrackHandle(ind)~=0)
            delete(hPlotTrackHandle(ind));
            a = find(ind ~= 0);
            for n = 1:length(a)
                trackedPos{a(n)} = [];
            end
            hPlotTrackHandle(ind) = 0;
        end
        
        newLight = step(tmLight, dt);
        if newLight ~= activeLight
            hLight(activeLight) = plot(ax(1), trafficLightPos(activeLight,1), trafficLightPos(activeLight,2),'o','MarkerSize',10,'MarkerFaceColor', 'w', 'MarkerEdgeColor','k');
            hLight(newLight) = plot(ax(1), trafficLightPos(newLight,1), trafficLightPos(newLight,2),'o','MarkerSize',10,'MarkerFaceColor', trafficLightColors(newLight), 'MarkerEdgeColor','k');
            activeLight = newLight;
        end
        
        for n=1:length(hStatGlobal)
            if(hStatGlobal(n) ~=0)
                delete(hStatGlobal(n));
            end
        end
        
        for lane = 1:4
            vehicleCount(lane) = vehicleCountArrived(lane) - stat(lane).pending;
        end
        
        hStatGlobal(1) = text(ax(1), -40,80, ['Frame #', num2str(gtick)]);
        hStatGlobal(2) = text(ax(1), -40,76, ['Total Vehicles: ', num2str(sum(vehicleCount)),' =[', num2str(vehicleCount(1)),', ',num2str(vehicleCount(2)),', ',num2str(vehicleCount(3)),', ',num2str(vehicleCount(4)),']']);
        hStatGlobal(3) = text(ax(1), -40,74, ['Total Objects:  ', num2str(sum(objectCount)),' =[', num2str(objectCount(1)),', ',num2str(objectCount(2)),', ',num2str(objectCount(3)),', ',num2str(objectCount(4)),']']);
        
        uicontrol(hRb);
        value = get(hRb, 'Value');
        while (value == 0)
            pause(10/1000);
            value = get(hRb, 'Value');
        end
        F(tick) = getframe(hFigure);
        pause(10/1000);
    end
    
    % Save one minute run
    for lane = 1:4
        vHist(lane).veh = [];
        vNum = 1;
        for ind=1:maxNumTracks
            marker = getVehicleMark(tmLane{lane},ind);
            if marker == 1
                hVeh = getVehicleHandle(tmLane{lane},ind);
                vHist(lane).veh{vNum} = hVeh;
                vNum = vNum+1;
                setVehicleMark(tmLane{lane},ind, 0);
            end
        end
        vHist(lane).num = vNum-1;
    end
    
    ind = (activeTracks == 3);
    tHist = tHist_current(ind);
    activeTracks(ind) = 0; %set back to free
 %{   
    tHist.track = [];
    tNum = 1;
    for ind = 1:maxNumTracks
        marker = getTrackerMark(hRadarTracking, ind);
        if marker == 1
            hTrack = getTrackerHandle(hRadarTracking, ind);
            tHist.track{tNum} = hTrack;
            tNum = tNum + 1;
            % clear marker
            setTrackerMark(hRadarTracking, ind, 0);
        end
        tHist.num = tNum-1;
    end
%}    
    disp('Saving...');
    matFilename = ['tm1500_', num2str(minuteTick,'%03d'),'.mat'];
    F1=F(1:150); F2=F(151:300); F3=F(301:450); F4=F(451:600); F5=F(601:750);  
    F6=F(751:900); F7=F(901:1050); F8=F(1051:1200); F9=F(1201:1350); F10=F(1351:1500);
    F = [];
    save(matFilename,'cHist','vHist', 'tHist','F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','-v7.3');
    F1=[]; F2=[]; F3=[]; F4=[]; F5=[]; F6=[]; F7=[]; F8=[]; F9=[]; F10=[];
end
disp(['Total Vehicles: ', num2str(sum(vehicleCount)),' =[', num2str(vehicleCount(1)),', ',num2str(vehicleCount(2)),', ',num2str(vehicleCount(3)),', ',num2str(vehicleCount(4)),']']);
disp(['Total Objects:  ', num2str(sum(objectCount)),' =[', num2str(objectCount(1)),', ',num2str(objectCount(2)),', ',num2str(objectCount(3)),', ',num2str(objectCount(4)),']']);
disp('Done');
