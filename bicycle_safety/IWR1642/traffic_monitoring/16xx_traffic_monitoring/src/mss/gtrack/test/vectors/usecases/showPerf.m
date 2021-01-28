tnum = 1;
vnum = ones(1,4);
vHistStruct = struct('veh', []);
vHistAll = repmat(vHistStruct,4, 1);
objectTrackCount = zeros(1,4);
filePath = 'C:\code\';
for n=1:10
    matFilename = [filePath,'tmRandomTestCFAR1NODW02UV100_', num2str(n,'%03d'),'.mat'];
    load(matFilename, 'tHist','vHist');
    for t=1:size(tHist.track, 2)
        tHistAll{tnum} = tHist.track{t};
        tnum = tnum +1;
    end
    for lane = 1:4
        for v=1:length(vHist(lane).veh)
            index = vHist(lane).veh{v}.id;
            vHistAll(lane).veh{index} = vHist(lane).veh{v};
        end
    end
end
numTracks = length(tHistAll);
trackingPrecisionPosition = cell(1,numTracks);
trackingPrecisionVelocity = cell(1,numTracks);
hFigure = figure('Name','Tracking Precision','NumberTitle','off', 'units','normalized','Position',[0 0 1 1]);
subplot(2,2,1);
xlabel('Range, m');
ylabel('X Error, m');
hold on;
subplot(2,2,2);
xlabel('Range, m');
ylabel('Y Error, m');
hold on;
subplot(2,2,3);
xlabel('Range, m');
ylabel('Vx Error, m/c');
hold on;
subplot(2,2,4);
xlabel('Range, m');
ylabel('Vy Error, m/c');
hold on;
hAxStat = gca;

hFigureError = figure('Name','Error Tracks','NumberTitle','off', 'units','normalized','Position',[0 0 1 1]);
hold on;
hAxError = gca;
hFigureStd = figure('Name','Std','NumberTitle','off');
hold on;
hAxStd = gca;

lbins = 15:1:70;
e{1}.stat=zeros(length(lbins),numTracks);
e{2}.stat=zeros(length(lbins),numTracks);
e{3}.stat=zeros(length(lbins),numTracks);
e{4}.stat=zeros(length(lbins),numTracks);
eMean =zeros(1,length(lbins));
eStd =zeros(1,length(lbins));
h = [];
legendStr = {};
hError = [];
legendErrorStr = {};
errorInsideThreshold = 1;
errorMaxThreshold = 3;

goodTrackCount = 0;
errorTrackCount = 0;
removedTrackCount = 0;

figure(hFigure);
for tid = 1:length(tHistAll)
    if(tid == 34)
        disp(tid);
    end
    hTrack = tHistAll{tid};
    if isempty(hTrack)
        break;
    end
    lane = floor(hTrack.vehicleId/1024);
    id = hTrack.vehicleId - lane*1024;

    if hTrack.tick < 20
%        hError = [hError, plot(hAxError, hTrack.sHatHistory(2:hTrack.tick,2),'.-')];
%        legendErrorStr = [legendErrorStr; ['TS: ',num2str(tid),': (',num2str(lane),'/',num2str(id), ')']];
%        legend(hError, legendErrorStr);
        disp(['tid=', num2str(tid), ' removed, too short']);
        removedTrackCount = removedTrackCount + 1;
        continue;
    end
    hVeh = vHistAll(lane).veh{id};
    %Sync the time stamps
    trackAllocationTime = hTrack.allocationTime;
    if(hTrack.tick > 400)
        tSkip = hTrack.tick - 400;
        tStart = mod(hTrack.tick, 400)+1;
        trackingTicksAll = [tStart:400, 1:tStart-1];
    else
        tSkip = 1;
        trackingTicksAll = 2:hTrack.tick;
    end
    vehicleEntryTime = hVeh.entryTime;
    shiftTS = trackAllocationTime - vehicleEntryTime + tSkip -1;

    trackingPosAll = hTrack.sHatHistory(trackingTicksAll,1:2);
    ind = (trackingPosAll(:,2) > 17); % Only interested in 17m+
    trackingTicks = trackingTicksAll(ind);
    if isempty(trackingTicks)
        disp(['tid=', num2str(tid), ' removed, no points >17m']);
        removedTrackCount = removedTrackCount +1;
        continue;
    end
    trackingPos = hTrack.sHatHistory(trackingTicks,1:2);    
    trackingVel = hTrack.sHatHistory(trackingTicks,3:4);
    if(hVeh.tick > 500)
        disp(['tid=', num2str(tid), ' removed, ground truth overwritten']);
        removedTrackCount = removedTrackCount +1;
        continue;
    end
    groundTruthTicks = 2:hVeh.tick;        
    groundTruthPosAll = hVeh.posxyHistory(groundTruthTicks,:) + [0,2];    
    groundTruthVelAll = hVeh.velxyHistory(groundTruthTicks,:);    
    ticks = 1:length(trackingTicks);
 
    minErrorBF = inf;
    for shift = 0:(length(groundTruthTicks)-length(trackingTicks))
        tickWindow = ticks+shift;
%        ind = abs(trackingVel(:,2)) > 0.1;
%        error = sum((groundTruthPosAll(tickWindow(ind),2)-trackingPos(ticks(ind),2)).^2);
        errorBF = sum((groundTruthPosAll(tickWindow,2)-trackingPos(ticks,2)).^2);
        if errorBF < minErrorBF
            minErrorBF = errorBF;
            shiftBF = shift;
        end
    end

%    groundTruthPosBF = groundTruthPosAll(ticks+shiftBF-1,:);
%    groundTruthVelBF = hVeh.velxyHistory(ticks+shiftBF-1,:);

    if(ticks(end)+shiftTS > size(groundTruthPosAll,1))
        disp(['tid=', num2str(tid), ' removed, tracked passed ground truth']);
        removedTrackCount = removedTrackCount +1;
        continue;
    else
        groundTruthPos = groundTruthPosAll(ticks+shiftTS,:);
        groundTruthVel = hVeh.velxyHistory(ticks+shiftTS,:);
        %ignore first 10
        errorVec = sqrt((groundTruthPos(5:end,1)-trackingPos(5:end,1)).^2 + 0.5*(groundTruthPos(5:end,2)-trackingPos(5:end,2)).^2);
        errorVecMean = mean(errorVec);
        numOutside = nnz(errorVec > errorInsideThreshold);
        if(max(errorVec) < errorMaxThreshold) && (numOutside/length(trackingPos(:,1)) < 0.1)
            goodTrackCount  = goodTrackCount + 1;
            resultBin = 1;
            resultString = 'Pass';
        else
            resultBin = 0;
            resultString = 'Fail';
            errorTrackCount = errorTrackCount +1;            
        end
    end
    if(shiftBF ~= shiftTS)
        disp(['Track: ',num2str(tid),', shift with BestFit is ', num2str(shiftBF), ', with TS is ', num2str(shiftTS)]); 
    end

%{    
    if(tid == 111)
        plot(hAxError, groundTruthPos(:,2),'.b-');
        plot(hAxError, trackingPos(:,2),'.r-');
    end
%}
    disp(['tid=',num2str(tid),'@',num2str(trackAllocationTime),': vid=',num2str(hTrack.vehicleId),'(',num2str(lane),'.',num2str(id),'), Error=', num2str(errorVecMean,3),' =>' resultString]);

    if(resultBin == 1)
        objectTrackCount(lane) = objectTrackCount(lane) +1;
        
        posErrorX = groundTruthPos(:,1)-trackingPos(:,1);
        posErrorY = groundTruthPos(:,2)-trackingPos(:,2);
        velErrorX = groundTruthVel(:,1)-trackingVel(:,1);
        velErrorY = -groundTruthVel(:,2)-trackingVel(:,2);
        if(nnz(abs(velErrorY) > 40))
            disp(tid);
        end
        
        trackingPrecisionPosition{tid}.ErrorX = posErrorX;
        trackingPrecisionPosition{tid}.ErrorY = posErrorY;
        trackingPrecisionVelocity{tid}.ErrorVx = velErrorX;
        trackingPrecisionVelocity{tid}.ErrorVy = velErrorY;
        
        [a,b,c] = histcounts(groundTruthPos(:,2), lbins);
        e{1}.stat(c+1,tid) = posErrorX;
        e{2}.stat(c+1,tid) = posErrorY;
        e{3}.stat(c+1,tid) = velErrorX;
        e{4}.stat(c+1,tid) = velErrorY;

        subplot(2,2,1);
        plot(groundTruthPos(:,2),posErrorX, '.b-');
        subplot(2,2,2);
%{        
        h = [h, plot(groundTruthPos(:,2),posErrorY)];
        legendStr = [legendStr; [num2str(tid),': (',num2str(lane),'/',num2str(id), ')']];
        legend(h, legendStr);
%}        
        plot(groundTruthPos(:,2),posErrorY);
        subplot(2,2,3);
        plot(groundTruthPos(:,2),velErrorX, '.r-');
        subplot(2,2,4);
        plot(groundTruthPos(:,2),velErrorY, '.r-');
    else        
        hError = [hError, plot(hAxError, groundTruthPos(:,2),'.-', 'color',hTrack.color)];
        legendErrorStr = [legendErrorStr; ['GT: ',num2str(tid),': (',num2str(lane),'/',num2str(id), ')']];
        hError = [hError, plot(hAxError, trackingPos(:,2),'.-','color',hTrack.color)];
        legendErrorStr = [legendErrorStr; ['TR: ',num2str(tid),': (',num2str(lane),'/',num2str(id), ')']];
        legend(hError, legendErrorStr);
        disp(['tid=', num2str(tid), ' removed, large error']);
    end
end
disp(['Vehicles tracked, Total =', num2str(sum(objectTrackCount)),',[',num2str(objectTrackCount(1)),',',num2str(objectTrackCount(2)),',',num2str(objectTrackCount(3)),',',num2str(objectTrackCount(4)),']']);
disp(['Total tracks =', num2str(length(tHistAll)), ',Good = ', num2str(goodTrackCount), ',Error = ', num2str(errorTrackCount), ',removed = ',num2str(removedTrackCount)]);

hStd=zeros(1,4);
eM=zeros(1,4);
eS=zeros(1,4);
descr = {'Xpos','Ypos','Vx','Vy'};
for splot=1:4
    subplot(2,2,splot);
    eStat = e{splot}.stat;
    for n=1:length(lbins)
        ind=eStat(n,:)~=0;
        if ~isempty(eStat(n,ind))
            eMean(n)=mean(eStat(n,ind));
            eStd(n)=std(eStat(n,ind));
        end
    end
    ind = eMean~=0;
    h1 = plot(lbins(ind), eMean(ind), 'LineWidth',2,'color','m');
    h2 = plot(lbins(ind), eStd(ind), 'LineWidth',2,'color','k');
    legend([h1, h2], 'Mean Error', 'Std');
    
    hStd(splot) = plot(hAxStd, lbins(ind), eStd(ind),'.-');
    
    eM(splot) = mean(eMean);
    eS(splot) = mean(eStd);
    disp([descr{splot},' Error, Mean =', num2str(eM(splot),2),', Std =',num2str(eS(splot),2)]);
end
legend(hAxStd, descr);
