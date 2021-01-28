colors='bgrcm';

GTRACK_PATH = 'C:\Radar\mmwave_sdk\ti\alg\gtrack';
GTRACK_MEX_PATH = [GTRACK_PATH,'\test\win\matlab\mex'];
GTRACK_MATLAB_PATH = [GTRACK_PATH,'\test\win\matlab\src'];
GTRACK_TESTVEC_PATH = [GTRACK_PATH,'\test\vectors\urv20'];

path(path,GTRACK_MEX_PATH);
path(path,GTRACK_MATLAB_PATH);
path(path,GTRACK_TESTVEC_PATH);

trackerRun = 'PCMex';
%trackerRun = 'PCMatlab';

maxNumTracks = 20;
maxNumPoints = 250;

%sensor parameters
sensor.rangeMax = 70;
sensor.rangeMin = 15;
sensor.azimuthFoV = 30*pi/180; %30 degree FOV in horizontal direction
sensor.azimuthTilt = asin(10/70);
sensor.azimuthMin = (pi/2 - (sensor.azimuthTilt + sensor.azimuthFoV/2));
sensor.azimuthMax = sensor.azimuthMin+sensor.azimuthFoV;
sensor.framePeriod = 40;
sensor.maxURadialVelocity = 20;

hStatGlobal = zeros(1,3);

tHistStruct = struct('tid', [], 'vid', [], 'allocationTime', 0, 'tick', 0, 'sHat', []);
tHist = repmat(tHistStruct, 1, 100);
hPlotCloudHandleAll = [];

%scene parameters
scene.numberOfLanes = 4;
scene.laneWidth = [3 3 3 3];
scene.laneXPos = [-3 0 3 6];
scene.stopLineDist = [21 20 20 20];
scene.stopExitDist = [15 15 15 15];
scene.trafficTimeouts = [20 2 5]; % [20 2 10] Green time, Yellow time, Red Time, (s)
scene.lambda = [1/10, 2/10, 3/10, 4/10]; % arrival rate per second, (1 vehicle in 10 seconds)
%scene.vehicleProb = [0.2 0.2 0.1 0.5; 0.3 0.6 0.1 0; 0.5 0.4 0.1 0; 0.9 0 0.1 0]; %Cars/Trucks/Motorcycles/Bicycles
scene.vehicleProb = [0.5 0 0 0.5; 1 0 0 0; 1 0 0 0; 1 0 0 0]; %Cars/Trucks/Motorcycles/Bicycles
angle=linspace(sensor.azimuthMin, sensor.azimuthMax, 128);

scrsz = get(groot,'ScreenSize');
hFigure = figure('Name','Traffic Monitoring','NumberTitle','off','Position',[0 0 1600 900]);

hRb = uicontrol(hFigure,'style','radio','string','Pause',...
                 'Position',[10 50 50 20],...
                 'Value',1,'UserData',1);
             
for n=1:4
    subplot(1,4,n);
    plot([0 sensor.rangeMax*cos(angle) 0],[0 sensor.rangeMax*sin(angle) 0], '-k'); hold on
    plot(sensor.rangeMin*cos(angle),sensor.rangeMin*sin(angle), '-k');    
    xlabel('X coordinate, m');
    ylabel('Y coordinate, m');
    if(n == 3)
        title('C/Mex library simulation');   
    end
    if(n == 4)
        title('Matlab library simulation');   
    end
    
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
trackingAx = ax(3);
gatingAx = ax(4);

trackingHistStruct = struct('tid', 0, 'allocationTime', 0, 'tick', 0, 'posIndex', 0, 'histIndex', 0, 'sHat', zeros(1000,6), 'ec', zeros(1000,9),'pos', zeros(100,2), 'hMeshU', [], 'hMeshG', [], 'hPlotTrack', [], 'hPlotCentroid', []);
trackingHist = repmat(trackingHistStruct, 1, maxNumTracks);

trackerConfig.stateVectorType = 1; % 0=>2D, 1=>2DA
trackerConfig.maxNumPoints = maxNumPoints;
trackerConfig.maxNumTracks = maxNumTracks;

trackerConfig.maxURadialVelocity = sensor.maxURadialVelocity;
trackerConfig.maxAcceleration = 20;
trackerConfig.deltaT = sensor.framePeriod/1000;
trackerConfig.initialRadialVelocity = -20;
trackerConfig.ax = [];

trackerConfig.advParams.allocation = [-1, 0, 3, 4, 2]; % any SNR, 1m/s doppler, a set of 3 Points with 4m in distance and 2m/s velocity separation
trackerConfig.advParams.gating = [16,12,8,0]; % 16 units of Volume, Limits: 12m in length, 8m in width and no limit it velocity
trackerConfig.advParams.thresholds = [3,3, 5,5,1]; % det2active, det2forget, act2free, stat2free, exit2free
trackerConfig.advParams.variations = [1/sqrt(12), 1/sqrt(12), 2]; % Height, Width, and Doppler std
trackerConfig.advParams.scenery = [-100,100,20,100]; % No walls, 10m lower exit zone, no upper zone. 

if(strcmp(trackerRun,'PCMex'))
    trackerConfig.verbose = 3;
    hRadarTrackingC = gtrack_moduleCreate_mex(trackerConfig);
elseif(strcmp(trackerRun,'PCMatlab'))
    trackerConfig.verbose = 3;
    hRadarTrackingM = trackModule('trackerConfig', trackerConfig);
end

activeTracks = zeros(1, trackerConfig.maxNumTracks);
for loop=1:10
    disp(['Loop #',num2str(loop)]);
    testVectorFile = ['pointCloud1500_', num2str(loop,'%03d'), '.dat'];
    tHist = repmat(tHistStruct, 1, 100);
    
    fCloud = fopen(testVectorFile, 'r');
    for tick=1:1500
        gtick = (loop -1)*1500 + tick;

        if(gtick == 190)
            disp(gtick)
        end

        numInputPoints = fscanf(fCloud, 'NUM=%d\n', 1);
        if(numInputPoints ~= 0)
            mRead = single(fscanf(fCloud, '%f,%f,%f,%f,%f,%f,%f\n',[7, numInputPoints]));
            pointCloud = mRead(1:4,:);
            
            % Roll radial velocity into +/- uRadVelocityMax space
            pointCloud(3,:) = mod(pointCloud(3,:)+sensor.maxURadialVelocity,2*sensor.maxURadialVelocity)-sensor.maxURadialVelocity;
            numOutputPoints = numInputPoints;
        else
            numOutputPoints = 0;
            pointCloud = single(zeros(3,0));
        end            

        switch trackerRun
            case 'PCMex'        
                [TID, S, EC, G] = gtrack_moduleStep_mex(hRadarTrackingC, pointCloud, [], numOutputPoints);
                    
            case 'PCMatlab'            
                [track, num] = step(hRadarTrackingM, pointCloud, 0, numOutputPoints);
                %translate from Matlab to C
                TID = zeros(1,num);
                S = zeros(6,num);
                EC = zeros(9,num);
                for n=1:num
                    TID(n) = track(n).trackerId -1; % In Matlab, TID are one-based
                    S(:,n) = track(n).S;
                    EC(:,n) = reshape(track(n).ec, 9,1);
                    G(n) = track(n).g;
                end            
        end
        
        if(ishandle(hPlotCloudHandleAll))
            delete(hPlotCloudHandleAll);
        end        
        if(numOutputPoints)                
            posAll = [pointCloud(1,:).*sin(pointCloud(2,:)); pointCloud(1,:).*cos(pointCloud(2,:))];            
            hPlotCloudHandleAll = plot(ax(4), posAll(1,:), posAll(2,:),'.k');
        end
        
        tNumC = length(TID);
        for n=1:tNumC
            tid = TID(n) + 1;
            g = G(n);
            centroid = computeH(S(:,n));
            ec = reshape(EC(:,n),3,3);
            if(nnz(ec)>1)
                [xU, yU, zU, vU] = gatePlot3(gatingAx, 1, centroid, ec);
                if isempty(trackingHist(tid).hMeshU)
                    trackingHist(tid).hMeshU = mesh(gatingAx, xU.*sin(yU),xU.*cos(yU), zU);
                    trackingHist(tid).hMeshU.EdgeColor = [0.5 0.5 0.5];
                    trackingHist(tid).hMeshU.FaceColor = 'none';
                else
                    set(trackingHist(tid).hMeshU, 'XData', xU.*sin(yU),'YData',xU.*cos(yU), 'ZData', zU);
                end
            end
            if(g ~= 0)
                [xG, yG, zG, vG] = gatePlot3(gatingAx, g, centroid, ec);
                if isempty(trackingHist(tid).hMeshG)
                    trackingHist(tid).hMeshG = mesh(gatingAx, xG.*sin(yG),xG.*cos(yG), zG);
                    trackingHist(tid).hMeshG.EdgeColor = colors(mod(tid,length(colors))+1);
                    trackingHist(tid).hMeshG.FaceColor = 'none';
                else
                    set(trackingHist(tid).hMeshG, 'XData', xG.*sin(yG),'YData',xG.*cos(yG), 'ZData', zG);
                end
            end
        end        

        for n=1:tNumC
            tid = TID(n) + 1;
            if(activeTracks(tid) == 0)
                activeTracks(tid) = 1;
                trackingHist(tid).tid = TID(n);
                trackingHist(tid).allocationTime = gtick;
                trackingHist(tid).tick = 1;
                trackingHist(tid).posIndex = 1;
                trackingHist(tid).histIndex = 1;
                trackingHist(tid).sHat(1,:) = S(:,n);
                trackingHist(tid).pos(1,:) = S(1:2,n);               
                trackingHist(tid).hPlotTrack = plot(trackingAx, S(1,n), S(2,n), '.-', 'color', colors(mod(tid,length(colors))+1));
                trackingHist(tid).hPlotCentroid = plot(trackingAx, S(1,n), S(2,n), 'o', 'color', colors(mod(tid,length(colors))+1));
            else
                activeTracks(tid) = 1;                
                trackingHist(tid).tick = trackingHist(tid).tick + 1;
                
                trackingHist(tid).histIndex = trackingHist(tid).histIndex + 1;
                if(trackingHist(tid).histIndex > 1000)
                    trackingHist(tid).histIndex = 1;
                end
                trackingHist(tid).sHat(trackingHist(tid).histIndex,:) = S(:,n);
                trackingHist(tid).ec(trackingHist(tid).histIndex,:) = EC(:,n);
                
                trackingHist(tid).posIndex = trackingHist(tid).posIndex + 1;
                if(trackingHist(tid).posIndex > 100)
                    trackingHist(tid).posIndex = 1;
                end
                trackingHist(tid).pos(trackingHist(tid).posIndex,:) = S(1:2,n);

                if(trackingHist(tid).tick > 100)
                    set(trackingHist(tid).hPlotTrack, 'XData', [trackingHist(tid).pos(trackingHist(tid).posIndex+1:end,1); trackingHist(tid).pos(1:trackingHist(tid).posIndex,1)], ...
                        'YData',[trackingHist(tid).pos(trackingHist(tid).posIndex+1:end,2); trackingHist(tid).pos(1:trackingHist(tid).posIndex,2)]);
                else
                    set(trackingHist(tid).hPlotTrack, 'XData', trackingHist(tid).pos(1:trackingHist(tid).posIndex,1), ...
                        'YData',trackingHist(tid).pos(1:trackingHist(tid).posIndex,2));
                end                
                set(trackingHist(tid).hPlotCentroid,'XData',S(1,n),'YData',S(2,n));                
            end
        end
        iDelete = find(activeTracks == 2);
        for n=1:length(iDelete)
            ind = iDelete(n);
            delete(trackingHist(ind).hPlotTrack);
            delete(trackingHist(ind).hPlotCentroid);
            delete(trackingHist(ind).hMeshU);
            delete(trackingHist(ind).hMeshG);
            trackingHist(ind).hMeshU = [];
            trackingHist(ind).hMeshG = [];
            activeTracks(ind) = 0;
        end
        
        iReady = (activeTracks == 1);
        activeTracks(iReady) = 2;
                
%{
        if(tNumM ~= tNum)
            disp(['Length Error: tNum(C)=', num2str(tNum), ', tNum(Matlab)=', num2str(tNumM)]);
        else
            if(tNum)
                for n=1:tNum
                    error = abs(S(:,n)-t(n).S);
                    ind = error > [0.1 0.1 0.1 0.1 1 1]';
                    if(ind)
                        disp(['Error at Tick ',num2str(gtick),',TID = ', num2str(TID(n))]); disp(S(:,n)); disp(t(n).S);
                    end
                end
            end
        end
%}
        
        hTextHandle = text(ax(1), -40,80, ['Frame #', num2str(gtick)]); 
        
        uicontrol(hRb);
        value = get(hRb, 'Value');
        while (value == 0)
            pause(10/1000);
            value = get(hRb, 'Value');
        end        
        pause(10/1000);
        delete(hTextHandle);
    end
    fclose(fCloud);
    save(['thistC_',num2str(loop,'%03d'),'.mat'], 'tHist');
end
disp('Done');


function [H] = computeH(s)
    posx = s(1); posy = s(2); velx = s(3); vely = s(4);
    range = sqrt(posx^2+posy^2);
    if posy == 0
        azimuth = pi/2;
    elseif posy > 0
        azimuth = atan(posx/posy);
    else
        azimuth = atan(posx/posy) + pi;
    end
    doppler = (posx*velx+posy*vely)/range;
    H = [range azimuth doppler]';
end

