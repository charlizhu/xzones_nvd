clear, clc, close all
colors = 'bgrcm';

%% SETUP GUI
[hSetup fhistPath files mode] = setup_replaygtrack();

%% Init & Define Variables
chirpParams = hSetup.params.dataPath;
frame_period_ms = hSetup.params.frameCfg.framePeriodicity;
scene = hSetup.scene;
scene.numberOfTargetBoxes = scene.numberOfLanes;
scene.azimuthTilt =  hSetup.angle*pi/180;
scene.lineX = [scene.leftLineX scene.leftLineX+cumsum(scene.laneWidth)];
scene.numberOfTargetBoxes = scene.numberOfLanes;
scene.targetBox = [scene.lineX(1:end-1); repmat(scene.stopLineY(1), size(scene.laneWidth)); scene.laneWidth; repmat(scene.startLineY-scene.stopLineY(1),size(scene.laneWidth))]';

if (mode ==1)
    trackerRun = 'Target';
elseif (mode ==2)
    trackerRun = 'PCMatlab';
end

pointCloudFilePath = fhistPath;
pointCloudFiles = files;
PointCloudDataSource = 'pointCloudCaptured';
captureVideo = hSetup.checkboxEnableVideo.Value;

maxNumTracks = 20;
maxNumPoints = 250;

hPlotCloudHandleAll = [];
hPlotCloudHandleOutRange = [];
hPlotCloudHandleClutter = [];
hPlotCloudHandleStatic = [];
hPlotCloudHandleDynamic =[];
hPlotPoints3D = [];

hTargetBoxHandle = zeros(scene.numberOfTargetBoxes,1);

clutterPoints = zeros(2,1);
activeTracks = zeros(1, maxNumTracks);
targetCountTotal = zeros(1,scene.numberOfLanes);
trackingHistStruct = struct('vid', 0, 'tid', 0, 'allocationTime', 0, 'tick', 0, 'targetCountState', 0, 'posIndex', 0, 'histIndex', 0, 'sHat', zeros(1000,6), 'ec', zeros(1000,9),'pos', zeros(100,2), 'hMeshU', [], 'hMeshG', [], 'hPlotAssociatedPoints', [], 'hPlotTrack', [], 'hPlotCentroid', []);
trackingHist = repmat(trackingHistStruct, 1, 20);
tHistAll = repmat(trackingHistStruct, 1, 100);
tHistInd = 1;
vehInd = 1;


if(strcmp(trackerRun,'PCMex') || strcmp(trackerRun,'PCMatlab'))
    GTRACK_PATH = '..\src\mss\gtrack';
    GTRACK_MEX_PATH = [GTRACK_PATH,'\test\win\matlab\mex'];
    GTRACK_MATLAB_PATH = [GTRACK_PATH,'\test\win\matlab\src'];
    path(path,GTRACK_MEX_PATH);
    path(path,GTRACK_MATLAB_PATH);
end
global flagExit
flagExit = 0;

%% Setup Figure

figHandle = figure('Name', 'Visualizer','tag','mainFigure')%,'CloseRequestFcn',@fig_closereq); 
clf(figHandle);
set(figHandle, 'WindowStyle','normal');
set(figHandle,'Name','Texas Instruments - People Counting','NumberTitle','off')    
set(figHandle,'currentchar',' ')         % set a dummy character
warning off MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame
jframe=get(figHandle,'javaframe');
set(figHandle, 'MenuBar', 'none');
set(figHandle, 'Color', [0 0 0]);
pause(0.00001);
set(jframe,'Maximized',1); 
pause(0.00001);

figureTitles = {'Statistics', 'Point Cloud', 'Gating and Association', 'Doppler Map', 'Control', 'Cumulative Cloud'};
figureGroup = [1, 2, 3, 4, 5, 4];
numFigures = size(figureTitles, 2);
hFigure = zeros(1,numFigures);

% Setup tab dimensions
hTabGroup(1) = uitabgroup(figHandle, 'Position', [0.0 0.5 0.19 0.5]);
hTabGroup(5) = uitabgroup(figHandle, 'Position', [0.0 0.0 0.19 0.5]);

hTabGroup(2) = uitabgroup(figHandle, 'Position', [0.19 0.0 0.27 1]);
hTabGroup(3) = uitabgroup(figHandle, 'Position', [0.46 0.0 0.27 1]);
hTabGroup(4) = uitabgroup(figHandle, 'Position', [0.73 0.0 0.27 1]);

% Populate tabs with plots
for iFig = 1:numFigures
    hFigure(iFig) = uitab(hTabGroup(figureGroup(iFig)), 'Title', figureTitles{iFig});
    if(strcmp(figureTitles{iFig},'Point Cloud') || strcmp(figureTitles{iFig},'Gating and Association'))
        %setup axes
        ax = axes('parent', hFigure(iFig));
        
        %copy axes from setup to figure
        fa = copyobj(hSetup.axes1.Children, ax);
        
        xlabel('X coordinate, m');
        ylabel('Y coordinate, m');
        axis equal;
        axis([hSetup.axes1.XLim hSetup.axes1.YLim]);
        grid on; hold on; grid minor;
              
        if(strcmp(figureTitles{iFig},'Point Cloud'))
            trackingAx = ax;
        end
        if(strcmp(figureTitles{iFig},'Gating and Association'))
            gatingAx = ax; 
        end
    end
    
    
    if(strcmp(figureTitles{iFig},'Statistics'))
        ax = axes('parent', hFigure(iFig));
        hStatGlobal(1) = text(0, 0.9, 'Frame # 0', 'FontSize',12);
        hStatGlobal(2) = text(0, 0.7, 'Detection Points: 0','FontSize',12);
        hStatGlobal(3) = text(0, 0.5, 'Target Count:  0','FontSize',12);
        hStatGlobal(4) = text(0, 0.3, 'In Box Count:  0','FontSize',12);
        axis off;
    end
    
    if(strcmp(figureTitles{iFig},'Control'))
        cFig = iFig;
        hRbPause = uicontrol(hFigure(cFig),'Style','radio','String','Pause','FontSize', 12,...
            'Units', 'normalized', 'Position',[0.1 0.1 0.3 0.1],'Value',0);
        hPbExit = uicontrol(hFigure(cFig),'Style', 'pushbutton', 'String', 'Exit','FontSize', 12,...
            'Units', 'normalized','Position', [0.6 0.1 0.3 0.1],'Callback', @exitPressFcn);
        setappdata(hPbExit, 'exitKeyPressed', 0);
    end
    
    if(strcmp(figureTitles{iFig},'Doppler Map'))
        ax = axes('parent', hFigure(iFig));
        line('Xdata',[0 0],'YData',hSetup.plotBound.Y, 'Color','k', 'LineStyle', '--', 'LineWidth', 0.5);
        line('Xdata',[-chirpParams.maxVelocity -chirpParams.maxVelocity],'YData',hSetup.plotBound.Y, 'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
        line('Xdata',[chirpParams.maxVelocity chirpParams.maxVelocity],'YData',hSetup.plotBound.Y, 'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
        axis([-20 chirpParams.maxVelocity+chirpParams.dopplerResolutionMps hSetup.plotBound.Y]);
        xlabel('Doppler, m/s');
        ylabel('Y coordinate, m');
        grid on;
        hold on;
        dopplerAx = gca;
    end
    if(strcmp(figureTitles{iFig},'Cumulative Cloud'))
        ax = axes('parent', hFigure(iFig));
        fa = copyobj(hSetup.axes1.Children, ax);       
        xlabel('X coordinate, m');
        ylabel('Y coordinate, m');
        axis equal;
        axis([hSetup.axes1.XLim hSetup.axes1.YLim]);
        grid on; hold on;
        cumCloudAx = gca;
    end 
end


%% Create tracking module

% If tracking is at PC, create tracking module
if(strcmp(trackerRun,'PCMex') || strcmp(trackerRun,'PCMatlab'))
    trackerConfig.verbose = 3; %debug info
    trackerConfig.stateVectorType = 1; % 0=>2D, 1=>2DA
    trackerConfig.maxNumPoints = maxNumPoints;
    trackerConfig.maxNumTracks = maxNumTracks;

    trackerConfig.maxRadialVelocity = chirpParams.maxVelocity; %sensor.maxRadialVelocity;
    trackerConfig.radialVelocityResolution = chirpParams.dopplerResolutionMps; %sensor.radialVelocityResolution;
    trackerConfig.maxAcceleration = str2double(hSetup.uitableKF.Data(1,2:3));
    trackerConfig.deltaT = hSetup.params.frameCfg.framePeriodicity/1000; %in sec
    trackerConfig.initialRadialVelocity = str2double(hSetup.uitableKF.Data(1,1));    
    trackerConfig.ax = [];
    
    trackerConfig.advParams.scenery =  [hSetup.box.numBoundary  reshape(hSetup.box.boundaryBox', 1, []) hSetup.box.numStatic reshape(hSetup.box.staticBox', 1, [])]; 
    trackerConfig.advParams.allocation = str2double(hSetup.uitableAllocation.Data(2:end)); % Total SNR, min Velocity, min Points, distance, velocity separation
    trackerConfig.advParams.gating = str2double(hSetup.uitableGating.Data(2:end)); % Volume, Length Width and Velocity Limits
    trackerConfig.advParams.thresholds = str2double(hSetup.uitableState.Data(2:end)); % det2act, det2free, act2free, stat2free, exit2free
    trackerConfig.advParams.variations = str2double(hSetup.uitableVariation.Data(2:end)); % Height, Width, Doppler std

    if(strcmp(trackerRun,'PCMex'))       
        hRadarTrackingC = gtrack_create_mex(trackerConfig);
    elseif(strcmp(trackerRun,'PCMatlab'))
        hRadarTrackingM = trackModule('trackerConfig', trackerConfig);
    end
    
end
pause(0.0001);
drawnow;

%%

point3D = [];
frameNum = 0;
fprintf('------------------\n');

previousPointCloud = zeros(4,0, 'single');


%% Loop through selected files and play back
if iscell(files)
    numFiles = size(files,2);
else
    numFiles =1;
    temp{1}= files;
    files=temp;
end

if captureVideo
    F(numFiles*1000) = struct('cdata',[],'colormap',[]);
end

for i=1:numFiles
    
    switch PointCloudDataSource
        case 'pointCloudCaptured'    
            fHistFileName = [fhistPath files{i}];
            if exist(fHistFileName, 'file' )
                load(fHistFileName,'fHist');
                disp(['Loading data from ', fHistFileName, ' ...']);
                numFrames = size(fHist,2);
            else
                break;
            end
    end

    %fSnapShot = repmat(videoFrame, numFrames,1);    
    
    for frameNum=1:numFrames

        if(get(hRbPause, 'Value') == 1)
            while(get(hRbPause, 'Value') == 1)
%                 if(getappdata(hPbStep, 'stepKeyPressed') == 1)
%                     setappdata(hPbStep, 'stepKeyPressed', 0);
%                     break;
%                 end
                pause(1);
            end
        end 
        
        if(strcmp(PointCloudDataSource,'pointCloudCaptured'))
            targetFrameNum = fHist(frameNum).targetFrameNum;
            numInputPoints = single(size(fHist(frameNum).pointCloud,2));
            numOutputPoints = numInputPoints;

            if(strcmp(trackerRun,'Target'))
                pointCloud = previousPointCloud;
                previousPointCloud = fHist(frameNum).pointCloud;
            else
                pointCloud = fHist(frameNum).pointCloud;
            end
                        
            if isempty(pointCloud)
                posAll = [];
                snrAll = [];
            else
                pointCloud(2,:)= pointCloud(2,:) + scene.azimuthTilt;
                posAll = [pointCloud(1,:).*sin(pointCloud(2,:)); pointCloud(1,:).*cos(pointCloud(2,:))];
                snrAll = pointCloud(4,:);
            end
            fHist(frameNum).numOutputPoints = numOutputPoints;
            fHist(frameNum).numInputPoints = numInputPoints;
        end
        
        % Delete previous points
        if(ishandle(hPlotCloudHandleAll))
            delete(hPlotCloudHandleAll);
        end
        if(size(posAll,2))
            % Plot all points
            hPlotCloudHandleAll = scatter(trackingAx, posAll(1,:), posAll(2,:),'.k','SizeData',snrAll*10);
            if(~isempty(dopplerAx))
                plot(dopplerAx, pointCloud(3,:), posAll(2,:), '.k');
            end
            if(~isempty(cumCloudAx))
                plot3(cumCloudAx, posAll(1,:), posAll(2,:), pointCloud(4,:),'.b');
            end
        end
        
        
        switch trackerRun
            case 'PCMex'        
                [TID, S, EC, Gain, mIndex] = gtrack_step_mex(hRadarTrackingC, pointCloud, [], numOutputPoints);
                
                numTargets = length(TID);
                
                fHist(frameNum).numTargets = numTargets;
                fHist(frameNum).targetList.numTargets = numTargets;
                fHist(frameNum).targetList.TID = TID;
                fHist(frameNum).targetList.S = S;
                fHist(frameNum).targetList.EC = EC;
                fHist(frameNum).targetList.G = Gain;
                fHist(frameNum).indexArray = mIndex;                    
                    
            case 'PCMatlab'            
                [track, num, mIndex] = step(hRadarTrackingM, pointCloud, 0, numOutputPoints);
                %translate from Matlab to C
                TID = zeros(1,num);
                S = zeros(6,num);
                EC = zeros(9,num);
                for n=1:num
                    TID(n) = track(n).trackerId -1; % In Matlab, TID are one-based
                    S(:,n) = track(n).S;
                    EC(:,n) = reshape(track(n).ec, 9,1);
                    Gain(n) = track(n).g;
                end            
                mIndex = mIndex -1;

            case 'Target'
                if isempty(fHist(frameNum).header)
                    continue;
                end
                if isempty(fHist(frameNum).targetList)
                    TID = [];
                    S = [];
                    EC = [];
                    Gain = [];
                    mIndex  = [];
                else
                    TID = fHist(frameNum).targetList.TID;
                    S = fHist(frameNum).targetList.S;
                    EC = fHist(frameNum).targetList.EC;
                    Gain = fHist(frameNum).targetList.G;
                    mIndex = fHist(frameNum).indexArray;
                end                
        end
                
        if nnz(isnan(S))
            disp('Error: S contains NaNs');
            continue;
        end
        if nnz(isnan(EC))
            disp('Error: EC contains NaNs');
            continue;
        end
        
        tNumC = length(TID);
        targetCountCurrent = tNumC;
        targetCountInBox = zeros(scene.numberOfTargetBoxes,1);
 
        if(size(mIndex,1)) 
            mIndex = mIndex + 1;
        end
        
        point3D = [posAll; pointCloud(3,:)];
        % Plot previous frame's 3D points       
        if(size(point3D,2))   
            if isempty(hPlotPoints3D)
                hPlotPoints3D = plot3(gatingAx, point3D(1,:), point3D(2,:), point3D(3,:),'.k');
            else
                set(hPlotPoints3D, 'XData', point3D(1,:),'YData', point3D(2,:), 'ZData', point3D(3,:));
            end
        end     
        
        for n=1:tNumC
            
            tid = TID(n)+1;            
            if(tid > maxNumTracks)
                disp('Error: TID is wrong');
                continue;
            end

            if( (size(mIndex,2) > 0) && (size(mIndex,2) == size(point3D,2)) )
                tColor = colors(mod(tid,length(colors))+1);
                ind = (mIndex == tid);
                if nnz(ind)
                    if isempty(trackingHist(tid).hPlotAssociatedPoints)
                        trackingHist(tid).hPlotAssociatedPoints = plot3(gatingAx, point3D(1,ind), point3D(2,ind), point3D(3,ind),'o', 'color', tColor);
                    else
                        if ishandle(trackingHist(tid).hPlotAssociatedPoints)
                            set(trackingHist(tid).hPlotAssociatedPoints, 'XData', point3D(1,ind),'YData', point3D(2,ind), 'ZData', point3D(3,ind));
                        end
                    end
                end
            end
               
            g = Gain(n);
            centroid = computeH(1, S(:,n));
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
                %            [gLim, a, b, c] = gateCreateLimLocal(8, ec, centroid, [4,1,0]);
                [xG, yG, zG, vG] = gatePlot3(gatingAx, g, centroid, ec);
                if isempty(trackingHist(tid).hMeshG)
                    trackingHist(tid).hMeshG = mesh(gatingAx, xG.*sin(yG),xG.*cos(yG), zG);
                    trackingHist(tid).hMeshG.EdgeColor = colors(mod(tid,length(colors))+1);
                    trackingHist(tid).hMeshG.FaceColor = 'none';
                else
                    set(trackingHist(tid).hMeshG, 'XData', xG.*sin(yG),'YData',xG.*cos(yG), 'ZData', zG);
                end
            end
            
            if(activeTracks(tid) == 0)
                activeTracks(tid) = 1;

                vehInd = vehInd + 1;
                trackingHist(tid).vid = vehInd;
                trackingHist(tid).tid = TID(n);
                trackingHist(tid).allocationTime = targetFrameNum;
                trackingHist(tid).tick = 1;
                trackingHist(tid).targetCountState = 0;
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
                if(~isempty(dopplerAx))
                    plot(dopplerAx, S(4,n), S(2,n), 'o', 'color', colors(mod(tid,length(colors))+1));
                end
                for nBoxes = 1:scene.numberOfTargetBoxes
                    if( (S(1,n) > scene.targetBox(nBoxes,1)) && (S(1,n) < (scene.targetBox(nBoxes,1) + scene.targetBox(nBoxes,3))) && ...
                            (S(2,n) > scene.targetBox(nBoxes,2)) && (S(2,n) < (scene.targetBox(nBoxes,2)+scene.targetBox(nBoxes,4))) )
                        targetCountInBox(nBoxes) = targetCountInBox(nBoxes) +1;
                    end                
                end
                % Target Counting Logic
                if(trackingHist(tid).targetCountState == 0)
                    if(S(2,n) > 25)
                        trackingHist(tid).targetCountState = 1;
                    end
                end
                if(trackingHist(tid).targetCountState == 1)
                    if(S(2,n) < 20)
                        if((S(1,n) > scene.lineX(1)) && (S(1,n) < scene.lineX(end)))
                            ind = find(scene.lineX < S(1,n), 1, 'last');
                            targetCountTotal(ind) = targetCountTotal(ind) +1;
                        end
                        trackingHist(tid).targetCountState = 2;
                    end
                end
            end
        end
        
        iDelete = find(activeTracks == 2);
        for n=1:length(iDelete)
            ind = iDelete(n);
            delete(trackingHist(ind).hPlotTrack);
            delete(trackingHist(ind).hPlotCentroid);
            delete(trackingHist(ind).hMeshU);
            delete(trackingHist(ind).hMeshG);
            delete(trackingHist(ind).hPlotAssociatedPoints);
            trackingHist(ind).hMeshU = [];
            trackingHist(ind).hMeshG = [];
            activeTracks(ind) = 0;
            tHistAll(tHistInd) = trackingHist(ind);
            tHistInd = tHistInd + 1;
        end
        
        iReady = (activeTracks == 1);
        activeTracks(iReady) = 2;
        
        string = [];
        string{1} = sprintf('Frame #: %d',targetFrameNum);
        string{2} = sprintf('Detection Points: %d',numOutputPoints);        
        
        string{3} = sprintf('Target Current [in Box]: %d [', targetCountCurrent);
        for nBoxes = 1:scene.numberOfTargetBoxes-1
            string{3} = [string{3}, sprintf('%d, ', targetCountInBox(nBoxes))];
        end
        string{3} = [string{3}, sprintf('%d]', targetCountInBox(scene.numberOfTargetBoxes))];
        
        string{4} = sprintf('Target Total [per Line]: %d, [', sum(targetCountTotal));
        for nLines = 1:scene.numberOfLanes-1
            string{4} = [string{4}, sprintf('%d, ', targetCountTotal(nLines))];
        end
        string{4} = [string{4}, sprintf('%d]', targetCountTotal(scene.numberOfLanes))];

        for s=1:length(hStatGlobal)
            set(hStatGlobal(s),'String',string{s});
        end
        
        for nBoxes = 1:scene.numberOfTargetBoxes
            if( (hTargetBoxHandle(nBoxes) ~= 0) && (ishandle(hTargetBoxHandle(nBoxes))))
                delete(hTargetBoxHandle(nBoxes));
            end
            if(targetCountInBox(nBoxes))
                hTargetBoxHandle(nBoxes) = rectangle(trackingAx, 'Position', scene.targetBox(nBoxes,:), 'EdgeColor','r', 'LineStyle', '-', 'LineWidth', 2);
            end
        end
        
        if (captureVideo)
            drawnow
        elseif (hSetup.popupmenuSpeed.Value == 1)
            pause(1/10);
            drawnow
        elseif (hSetup.popupmenuSpeed.Value == 2)
            pause(1/100);
        else
            pause(1/10000);
        end
        
        if captureVideo
            F(frameNum+(i-1)*1000)= getframe(figHandle);
        end
        
        if(flagExit)
            if captureVideo
                disp('Saving video frames before quitting...')
                save([fhistPath '\tmFrames.mat'], 'F','-v7.3');
                disp('Save complete. Quitting...')
            end
                delete(figHandle)
                delete(hSetup.figure1)
            break;
        end
    end

end

if captureVideo
    save([fhistPath '\tmFrames.mat'], 'F', 'frame_period_ms','-v7.3');
end

tHist = tHistAll(1:tHistInd);
tReady = find(activeTracks == 1);
for ind = 1:length(tReady)
    tid = tReady(n);
    tHist(tHistInd) = trackingHist(tid);    
    tHistInd = tHistInd + 1;
end
save('tHist.mat','tHist');
disp('Done');


%% UTILS
function fig_closereq(src,callbackdata)
global flagExit
% Close request function 
% to display a question dialog box 
   selection = questdlg('Close This Figure?',...
      'Close Request Function',...
      'Yes','No','Yes'); 
   switch selection, 
      case 'Yes',
         flagExit = 1;
         disp('Sending Exit command')        
      case 'No'
      return 
   end
end

%Display Chirp parameters in table on screen
function h = displayChirpParams(Params, Position)

    dat =  {'Start Frequency (Ghz)', Params.profileCfg.startFreq;...
            'Slope (MHz/us)', Params.profileCfg.freqSlopeConst;...   
            'Samples per chirp', Params.profileCfg.numAdcSamples;...
            'Chirps per frame',  Params.dataPath.numChirpsPerFrame;...
            'Frame duration (ms)',  Params.frameCfg.framePeriodicity;...
            'Sampling rate (Msps)', Params.profileCfg.digOutSampleRate / 1000;...
            'Bandwidth (GHz)', Params.profileCfg.freqSlopeConst * Params.profileCfg.numAdcSamples /...
                               Params.profileCfg.digOutSampleRate;...
            'Range resolution (m)', Params.dataPath.rangeResolutionMeters;...
            'Velocity resolution (m/s)', Params.dataPath.dopplerResolutionMps;...
            'Number of Rx (MIMO)', Params.dataPath.numRxAnt; ...
            'Number of Tx (MIMO)', Params.dataPath.numTxAnt;};
    columnname =   {'Chirp Parameter (Units)      ', 'Value'};
    columnformat = {'char', 'numeric'};
    
    h = uitable('Units','normalized', ...
            'Position', Position, ...
            'Data', dat,... 
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnWidth', 'auto',...
            'RowName',[]);
end

function [P] = parseCfg(cliCfg)
    P=[];
    for k=1:length(cliCfg)
        C = strsplit(cliCfg{k});
        if strcmp(C{1},'channelCfg')
            P.channelCfg.txChannelEn = str2double(C{3});
            P.dataPath.numTxAzimAnt = bitand(bitshift(P.channelCfg.txChannelEn,0),1) +...
                                      bitand(bitshift(P.channelCfg.txChannelEn,-1),1);
            P.dataPath.numTxElevAnt = 0;
            P.channelCfg.rxChannelEn = str2double(C{2});
            P.dataPath.numRxAnt = bitand(bitshift(P.channelCfg.rxChannelEn,0),1) +...
                                  bitand(bitshift(P.channelCfg.rxChannelEn,-1),1) +...
                                  bitand(bitshift(P.channelCfg.rxChannelEn,-2),1) +...
                                  bitand(bitshift(P.channelCfg.rxChannelEn,-3),1);
            P.dataPath.numTxAnt = P.dataPath.numTxElevAnt + P.dataPath.numTxAzimAnt;
                                
        elseif strcmp(C{1},'dataFmt')
        elseif strcmp(C{1},'profileCfg')
            P.profileCfg.startFreq = str2double(C{3});
            P.profileCfg.idleTime =  str2double(C{4});
            P.profileCfg.rampEndTime = str2double(C{6});
            P.profileCfg.freqSlopeConst = str2double(C{9});
            P.profileCfg.numAdcSamples = str2double(C{11});
            P.profileCfg.digOutSampleRate = str2double(C{12}); %uints: ksps
        elseif strcmp(C{1},'chirpCfg')
        elseif strcmp(C{1},'frameCfg')
            P.frameCfg.chirpStartIdx = str2double(C{2});
            P.frameCfg.chirpEndIdx = str2double(C{3});
            P.frameCfg.numLoops = str2double(C{4});
            P.frameCfg.numFrames = str2double(C{5});
            P.frameCfg.framePeriodicity = str2double(C{6});
        elseif strcmp(C{1},'guiMonitor')
            P.guiMonitor.detectedObjects = str2double(C{2});
            P.guiMonitor.logMagRange = str2double(C{3});
            P.guiMonitor.rangeAzimuthHeatMap = str2double(C{4});
            P.guiMonitor.rangeDopplerHeatMap = str2double(C{5});
        end
    end
    P.dataPath.numChirpsPerFrame = (P.frameCfg.chirpEndIdx -...
                                            P.frameCfg.chirpStartIdx + 1) *...
                                            P.frameCfg.numLoops;
    P.dataPath.numDopplerBins = P.dataPath.numChirpsPerFrame / P.dataPath.numTxAnt;
    P.dataPath.numRangeBins = pow2roundup(P.profileCfg.numAdcSamples);
    P.dataPath.rangeResolutionMeters = 3e8 * P.profileCfg.digOutSampleRate * 1e3 /...
                     (2 * P.profileCfg.freqSlopeConst * 1e12 * P.profileCfg.numAdcSamples);
    P.dataPath.rangeIdxToMeters = 3e8 * P.profileCfg.digOutSampleRate * 1e3 /...
                     (2 * P.profileCfg.freqSlopeConst * 1e12 * P.dataPath.numRangeBins);
    P.dataPath.dopplerResolutionMps = 3e8 / (2*P.profileCfg.startFreq*1e9 *...
                                        (P.profileCfg.idleTime + P.profileCfg.rampEndTime) *...
                                        1e-6 * P.dataPath.numDopplerBins * P.dataPath.numTxAnt);
end

function stepPressFcn(hObject, ~)
    setappdata(hObject, 'stepKeyPressed', 1);
end
function exitPressFcn(hObject, ~)
    setappdata(hObject, 'exitKeyPressed', 1);
    global flagExit
    flagExit =1;
end

function clearCloudPressFcn(hObject, ~)
    sensor = getappdata(hObject, 'sensor');
    scene = getappdata(hObject, 'scene');
    cla;
    setCumCloudFigure(sensor, scene);
end

function clearDopplerPressFcn(hObject, ~)
    sensor = getappdata(hObject, 'sensor');
    cla;
    setDopplerMapFigure(sensor);
end

function setDopplerMapFigure(sensor)
    line('Xdata',[0 0],'YData',[0 80], 'Color','k', 'LineStyle', '--', 'LineWidth', 0.5);
    line('Xdata',[-sensor.maxRadialVelocity -sensor.maxRadialVelocity],'YData',[0 80], 'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    line('Xdata',[sensor.maxRadialVelocity sensor.maxRadialVelocity],'YData',[0 80], 'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    axis([-20 20 0 80]);
    xlabel('Doppler, m/s');
    ylabel('Y coordinate, m');
    grid on;
    hold on;
end

function setCumCloudFigure(sensor, scene)
    plot(sensor.rangeMin*sin(sensor.angles+scene.azimuthTilt), sensor.rangeMin*cos(sensor.angles+scene.azimuthTilt), '-k');  hold on;
    plot([0 sensor.rangeMax*sin(sensor.angles+scene.azimuthTilt) 0],[0 sensor.rangeMax*cos(sensor.angles+scene.azimuthTilt) 0], '-k');
    xlabel('X coordinate, m');
    ylabel('Y coordinate, m');
    axis(scene.maxPos);
    grid on;
end

function config = readCfg(filename)
    config = cell(1,100);
    fid = fopen(filename, 'r');
    if fid == -1
        fprintf('File %s not found!\n', filename);
        return;
    else
        fprintf('Opening configuration file %s ...\n', filename);
    end
    tline = fgetl(fid);
    k=1;
    while ischar(tline)
        config{k} = tline;
        tline = fgetl(fid);
        k = k + 1;
    end
    config = config(1:k-1);
    fclose(fid);
end

function [H] = computeH(~, s)
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

function [XX, YY, ZZ, v] = gatePlot3(~, G, C, A)
    %Extract the ellipsoid's axes lengths (a,b,c) and the rotation matrix (V) using singular value decomposition:
    [~,D,V] = svd(A/G);

    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));
    c = 1/sqrt(D(3,3));
    v = 4*pi*a*b*c/3;

    % generate ellipsoid at 0 origin
    [X,Y,Z] = ellipsoid(0,0,0,a,b,c);
    XX = zeros(size(X));
    YY = zeros(size(X));
    ZZ = zeros(size(X));
    for k = 1:length(X)
        for j = 1:length(X)
            point = [X(k,j) Y(k,j) Z(k,j)]';
            P = V * point;
            XX(k,j) = P(1)+C(1);
            YY(k,j) = P(2)+C(2);
            ZZ(k,j) = P(3)+C(3);
        end
    end
end

function [y] = pow2roundup (x)
    y = 1;
    while x > y
        y = y * 2;
    end
end
