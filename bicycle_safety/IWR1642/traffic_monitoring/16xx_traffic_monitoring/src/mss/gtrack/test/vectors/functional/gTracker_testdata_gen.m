%{
 *   @file  gTracker_testdata_gen.m
 *
 *   @brief
 *      %generate test data for gtracker
 *
 *  
 *  NOTE:
 *      (C) Copyright 2016 Texas Instruments, Inc.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *    Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 *    Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 *    Neither the name of Texas Instruments Incorporated nor the names of
 *    its contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
%}


clear
x = 8647;
rand('state',x);
randn('state',x);

filePtr = fopen('gtracker_testdata.dat','w');
genVarFlag = 0;   % if needed to write measurement variance in the test vector
numTests = 6;  % number of test
fprintf(filePtr,'%u\n',numTests);
% 1: 1 tracker, eclipse
% 2: 2 trackers
% 3: 1 linear tracker with accelaration
% 4: 1 tracker with different number of observations per frame
% 5: 2 trackers merge and then diverge
% 6: 1 tracker with deaccelaration

numTrackers = [1; 2;1;1;2;1];  
T = 0.4; 
numData = 20;   % per frame per tracker
numFrame = 40;  % number of frames
axy = [0; 0];  %accelation xy
F=[1 0 T 0 0.5*T^2 0; 0 1 0 T 0 0.5*T^2; 0 0 1 0 T 0; 0 0 0 1 0 T; 0 0 0 0 1 0; 0 0 0 0 0 1];


fprintf(filePtr,'%u\n',numFrame);

%generate different circle
minTheta = 0;
maxTheta = pi*2;
theta = minTheta:(maxTheta-minTheta)/numFrame:maxTheta;
linearFlag = 0;
for ii=1:numTests
     linearFlag  = 0;
     cloudeP=[];
     xloc=[];
     yloc=[];
     vx=[];
     vy=[];
     tempv=[];
    for kk = 1:numTrackers(ii)
        if kk==1   
            if numTrackers(ii)==1
                if ii==3 %linear
                    xy =[10;5];
                    vxy=[2;3.5];
                    deltaxy(1,:) = [0.8 0.8];  % the xy range
                    axy = [0.2;0.1];
                    xloc(kk,1) = xy(1);
                    yloc(kk,1) = xy(2);
                    vx(kk,1) = vxy(1);
                    vy(kk,1) = vxy(2);
                    linearFlag = 1;
                elseif ii==6  % linear
                    xy =[10;10];
                    vxy=[5;7];
                    deltaxy(1,:) = [1 1];  % the xy range
                    axy = [-0.05;-0.08];
                    xloc(kk,1) = xy(1);
                    yloc(kk,1) = xy(2);
                    vx(kk,1) = vxy(1);
                    vy(kk,1) = vxy(2);
                    linearFlag = 1;   
                elseif ii==4
                    numData = 30;
                    %xloc(kk,:) = [10:0.3:30];
                    xloc(kk,:) = [10:0.6:40];
                    a = 40;
                    %yloc(kk,:) =8*sin(2*pi*0.1*xloc(kk,:))+15;
                    yloc(kk,:)=8*sin(2*pi*0.05*xloc(kk,:))+15;
                    vx(kk,:) = diff(xloc(kk,:));
                    vy(kk,:) = diff(yloc(kk,:));
                    deltaxy(1,:) = [0.5 0.5];  % the xy range
                   
                elseif ii==1%circle
                    center = [18;18];
                    radius = [12;10];
                    deltaxy(1,:) = [0.8 0.8];  % the xy range
                    xloc(kk,:) = radius(1)*cos(theta)+center(1);
                    yloc(kk,:) = radius(2)*sin(theta)+center(2);
                    vx(kk,:) = diff(xloc(kk,:));
                    vy(kk,:) = diff(yloc(kk,:));
        
                end
            else
                if ii==2 %the second one
                center = [28;28];
                radius = [12;12];
                deltaxy(1,:) = [0.8 0.8];  % the xy range
                xloc(kk,:) = radius(1)*cos(theta)+center(1);
                yloc(kk,:) = radius(2)*sin(theta)+center(2);
                
                vx(kk,:) = diff(xloc(kk,:));
                vy(kk,:) = diff(yloc(kk,:));
                elseif ii==5
                    % 2 sin wave 
                    xloc(kk,:) = [20:0.6:50];
                    a = 40;
                    %yloc(kk,:) =8*sin(2*pi*0.1*xloc(kk,:))+15;
                    yloc(kk,:)=6*sin(2*pi*0.05*xloc(kk,:))+20;
                    vx(kk,:) = diff(xloc(kk,:));
                    vy(kk,:) = diff(yloc(kk,:));
                    deltaxy(1,:) = [0.6 0.6];  % the xy range
                end
        
            end
        
        else       %triangle
            if (ii==2) % the second one
            t=[-numFrame:2:numFrame];
            w=2*numFrame;
            u=tripuls(t,w);
            xloc(kk,:) = (t+numFrame)/4;
            yloc(kk,:) = u*20+ 5;
            
            deltaxy(2,:)= [0.3 0.3];  % the xy range
            
            vx(kk,:) = diff(xloc(kk,:));
            vy(kk,:) = diff(yloc(kk,:));
            else
               tx = [20:0.6:100];
               ty = 6*sin(2*pi*0.057*tx+pi*3.5/4)+20;
               
                 xloc(kk,:) = tx(1:length(xloc(1,:)));
                 yloc(kk,:)= ty(1:length(xloc(1,:)));
                 vx(kk,:) = diff(xloc(kk,:));
                 vy(kk,:) = diff(yloc(kk,:));
                 deltaxy(1,:) = [0.6 0.6];  % the xy range
            end
            
        end
        
       
        
       
        %generate xy location
        cloudeP(1,(kk-1)*numData+1:kk*numData) = xloc(kk,1)-deltaxy(kk,1)+2*deltaxy(kk,1)*rand(1,numData);    % x-location
        cloudeP(2,(kk-1)*numData+1:kk*numData) = yloc(kk,1)-deltaxy(kk,2)+2*deltaxy(kk,2)*rand(1,numData);    % y-location
        if (linearFlag ==1)
            cloudeP(3,(kk-1)*numData+1:kk*numData) = vx(kk,1)*ones(1,numData);   %x-velocity
            cloudeP(4,(kk-1)*numData+1:kk*numData) = vy(kk,1)*ones(1,numData);   %y-velocity
        else
        cloudeP(3,(kk-1)*numData+1:kk*numData) = vx(kk,1)*ones(1,numData)/T;   %x-velocity
        cloudeP(4,(kk-1)*numData+1:kk*numData) = vy(kk,1)*ones(1,numData)/T;   %y-velocity
        end
        cloudeP(5,(kk-1)*numData+1:kk*numData) = axy(1)*ones(1,numData);   %x-a
        cloudeP(6,(kk-1)*numData+1:kk*numData) = axy(2)*ones(1,numData);   %y-a
        
    end
    
    %generate the data
    cloudData=[];
    measureData=[];
    range=[];
    aoa=[];
    vel=[];
    varr=[];
    vara=[];
    varv=[];
    tempv=[];
    velstd=[];
    xstd=[];
    ystd=[];
    figure;
   
    for kk=1:numFrame
        
        if (ii==4)
            numdataperframe = randi([10,30],1,1); 
        else
            numdataperframe = numData*numTrackers(ii);
        end
        %numdataperframe
        fprintf(filePtr,'NUM=%d\n',numdataperframe);
         
         for jj = 1:numData*numTrackers(ii)
             range(jj) = sqrt(cloudeP(1,jj)^2+cloudeP(2,jj)^2);
             aoa(jj) = atan(cloudeP(1,jj)/cloudeP(2,jj));
             vel(jj) = (cloudeP(1,jj)*cloudeP(3,jj)+cloudeP(2,jj)*cloudeP(4,jj))/range(jj);
         end
         
        
         varrtemp=[];
        %snapnow;
         varatemp=[];
         varvtemp=[];
         
         for jj=1:numTrackers(ii)
             varr(jj) = 8*sum((range((jj-1)*numData+1:jj*numData)-mean(range((jj-1)*numData+1:jj*numData))).^2)/(numData-1);
             vara(jj) = 8*sum((aoa((jj-1)*numData+1:jj*numData)-mean(aoa((jj-1)*numData+1:jj*numData))).^2)/(numData-1);
             varv(jj) = 8*sum((vel((jj-1)*numData+1:jj*numData)-mean(vel((jj-1)*numData+1:jj*numData))).^2)/(numData-1);
          %dataTTempX = cloudeP(1,(jj-1)*numData+1:jj*numData);
          %dataTTempY = cloudeP(2,(jj-1)*numData+1:jj*numData);
          %centerXY = [mean(dataTTempX) mean(dataTTempY)];
          
             varrtemp = [varrtemp  varr(jj)*ones(1,numData)];
             varatemp= [varatemp vara(jj)*ones(1,numData)];
             varvtemp = [varvtemp  varv(jj)*ones(1,numData)];
             
             %find out the standard deviation
            if (1)
                velstd((kk-1)*numTrackers(ii)+jj) = max(vel((jj-1)*numData+1:jj*numData))-min(vel((jj-1)*numData+1:jj*numData));
                xstd((kk-1)*numTrackers(ii)+jj) = max(cloudeP(1,(jj-1)*numData+1:jj*numData))-min(cloudeP(1,(jj-1)*numData+1:jj*numData));
                ystd((kk-1)*numTrackers(ii)+jj) = max(cloudeP(1,(jj-1)*numData+1:jj*numData))-min(cloudeP(1,(jj-1)*numData+1:jj*numData));
            end

        
         
         end
         %permutate the data
         if (ii==4)
             pIndex=[1:numTrackers(ii) * numData];
         else
             pIndex = randperm(numTrackers(ii) * numData);
         end
         
        for jj = 1:numdataperframe
            tIdx = floor((pIndex(jj)-1)/numData)+1;
            if (genVarFlag==1)
                fprintf(filePtr,'%f,%f,%f,%f,%f,%f\n',range(pIndex(jj)),aoa(pIndex(jj)),vel(pIndex(jj)),varr(tIdx),vara(tIdx),varv(tIdx)); 
            else
                fprintf(filePtr,'%f,%f,%f\n',range(pIndex(jj)),aoa(pIndex(jj)),vel(pIndex(jj)));
            end
        
        end
    
        
        %save cloude
        cloudData(kk,:,:) = cloudeP;
        measureData(kk,1,:) = range;
        measureData(kk,2,:) = aoa;
        measureData(kk,3,:) = vel;
        measureData(kk,4,:) = varrtemp;
        measureData(kk,5,:) = varatemp;
        measureData(kk,6,:) = varvtemp;
        tempv=[tempv vel];
         
        
        
        if (kk==1)
             plot(squeeze(cloudData(1,1,1:numdataperframe)),squeeze(cloudData(1,2,1:numdataperframe)),'b.')
             hold on
            
        else
            if (mod(kk,2)==0)
                plot(squeeze(cloudData(kk,1,1:numdataperframe)),squeeze(cloudData(kk,2,1:numdataperframe)),'r.')
            else
                 plot(squeeze(cloudData(kk,1,1:numdataperframe)),squeeze(cloudData(kk,2,1:numdataperframe)),'b.')
            end
        end
         %add cloud mean as reference
        fprintf(filePtr,'%u\n',numTrackers(ii));
        for jj=1:numTrackers(ii)
            if (ii==4)
                cloudeRef = mean(cloudeP(:,(jj-1)*numData+1:(jj-1)*numData+numdataperframe),2);
                 vrRef = mean( measureData(kk,3,(jj-1)*numData+1:(jj-1)*numData+numdataperframe));
            else
                cloudeRef = mean(cloudeP(:,(jj-1)*numData+1:jj*numData),2);
                vrRef = mean( measureData(kk,3,(jj-1)*numData+1:jj*numData));
            end
            fprintf(filePtr,'%f,%f,%f,%f\n',cloudeRef(1),cloudeRef(2),cloudeRef(3),cloudeRef(4));
            plot(cloudeRef(1),cloudeRef(2),'k*');
           %vrRef = (cloudeRef(1)*cloudeRef(3)+cloudeRef(2)*cloudeRef(4))/sqrt(cloudeRef(1)*cloudeRef(1)+cloudeRef(2)*cloudeRef(2));
           %fprintf(filePtr,'%f,%f,%f\n',cloudeRef(1),cloudeRef(2),vrRef);
        end
        
         
        
        %update to next frame
        
        for jj=1:numTrackers(ii)
          
            temp = F*mean(cloudeP(:,(jj-1)*numData+1:jj*numData),2);
            
            cloudeP(1,(jj-1)*numData+1:jj*numData) = temp(1)-deltaxy(jj,1)+2*deltaxy(jj,1)*rand(1,numData);    % x-location
            cloudeP(2,(jj-1)*numData+1:jj*numData) = temp(2)-deltaxy(jj,2)+2*deltaxy(jj,2)*rand(1,numData);    % y-location 
        
            cloudeP(5,(jj-1)*numData+1:jj*numData) = axy(1)*ones(1,numData);   %x-a
            cloudeP(6,(jj-1)*numData+1:jj*numData) = axy(2)*ones(1,numData);   %y-a
            
            if linearFlag==1
                cloudeP(3,(jj-1)*numData+1:jj*numData) = temp(3)*ones(1,numData);   %x-velocity
                cloudeP(4,(jj-1)*numData+1:jj*numData) = temp(4)*ones(1,numData); 
            else
                
                 if (kk<numFrame)
                    cloudeP(3,(jj-1)*numData+1:jj*numData) = vx(jj,kk)*ones(1,numData)/T;   %x-velocity
                    cloudeP(4,(jj-1)*numData+1:jj*numData) = vy(jj,kk)*ones(1,numData)/T;   %y-velocity
                 end
            end
        end
       
        
        
    end
    
    axis equal
    if (0)
    figure;
    plot(tempv,'.');
    ylabel('veloctiy measurement');
    end
    if (0) %find the maximum acceleration
        maxa=[];
        temppv = reshape(tempv,numData*numTrackers(ii),numFrame);
        for mmm=1:numTrackers(ii)
            for jj=1:numFrame-1
               % tempa(:,jj) = temppv((mmm-1)*numData+1:mmm*numData,jj+1)-temppv((mmm-1)*numData+1:mmm*numData,jj);
               tempa(jj) = max(temppv((mmm-1)*numData+1:mmm*numData,jj+1))-min(temppv((mmm-1)*numData+1:mmm*numData,jj));
            end
            %maxa(mmm) = max(max(abs(tempa)))/T;
            maxa(mmm) = max(abs(tempa))/T;
        end
        maxa
    end
    
    if (0)
        figure;plot(velstd,'.');title('velocity max-min')
        figure;plot(xstd,'.');title('x max-min');
        figure;plot(ystd,'.');title('y max-min');
    end
end

 fclose(filePtr);
 filePtr = fopen('gtracker_benchmarkdata.dat','w');
% for benchmark data
% generate data for step function benchmark
% number of data per frame is 300, number of trackers is 20
for jj = 1:20  
    numdataPerObj = ceil(50/jj);
    if numdataPerObj<4  %in code allocation threshold is set to 3
        numdataPerObj = 4;
    end
    numdata = numdataPerObj *jj * 6;   % test 6 level of number of data
    if (jj<17)
        modIdx = 0:1:15;
        modIdxStep = floor(16/jj);
    else
        modIdx = 0:1:31;
        modIdxStep = floor(32/jj);
    end
    
    index=[];
    index = zeros(1,numdata); 
    for ii=0:jj-1
         index(ii*(numdata/jj)+1:(ii+1)*(numdata/jj))=modIdx(ii*modIdxStep+1);
    end
    
    if (jj<17)
        symbolTemp = qammod(index,16) * 15;
        symbolTempN = qammod(index,16,4/180*pi)*15;
    else
        symbolTemp = qammod(index,32) * 20;
        symbolTempN = qammod(index,32,4/180*pi)*20;
    end
    
    
    symbol = awgn(symbolTemp,15);
    
    if (0) %(jj==5)
        figure;
        plot(real(symbol),imag(symbol),'r.')
        
    end
    
    xsymbol = real(symbol);
    ysymbol = imag(symbol);
    symbolN = awgn(symbolTempN,15);
    xsymbolN = real(symbolN);
    ysymbolN = imag(symbolN);
    range=[];
    aoa=[];
    vel=[];
    for ii=1:numdata
       range(ii) = sqrt(xsymbol(ii)^2+ysymbol(ii)^2);
       if (ysymbol(ii)>0)
          aoa(ii) = atan(xsymbol(ii)/ysymbol(ii));
       else
          aoa(ii) = atan(xsymbol(ii)/ysymbol(ii))+pi;
       end
       
       velx = xsymbolN(ii)-xsymbol(ii);
       vely = ysymbolN(ii)-ysymbol(ii);
       vel(ii) = (xsymbol(ii)*velx+ysymbol(ii)*vely)/range(ii);

    end
    
    rangeT=[];
    aoaT=[];
    
    for ii=1:numdata
       rangeT(ii) = sqrt(xsymbolN(ii)^2+ysymbolN(ii)^2);
       if (ysymbol(ii)>0)
          aoaT(ii) = atan(xsymbolN(ii)/ysymbolN(ii));
       else
          aoaT(ii) = atan(xsymbolN(ii)/ysymbolN(ii))+pi;
       end
       
      

    end
    
    
    varr = [];
    vara = [];
    varv = [];
    varr = 4*sum((range-mean(range)).^2)/(numdata-1);
    vara = 4*sum((aoa-mean(aoa)).^2)/(numdata-1);
    varv = 4*sum((vel-mean(vel)).^2)/(numdata-1);
    %numdata
    fprintf(filePtr,'NUM=%d\n',numdata);
    for ii=1:6  %every time write 6 groups of numdataPerObj data
        for nn=1:jj  % the tracker index
            temp = (nn-1)* numdataPerObj * 6 +1 + (ii-1)*numdataPerObj;
            for mm=0:numdataPerObj-1
               % temp+mm
               if (genVarFlag==1)
                   fprintf(filePtr,'%f,%f,%f,%f,%f,%f\n',range(temp+mm),aoa(temp+mm),vel(temp+mm),varr,vara,varv);
               else
                   fprintf(filePtr,'%f,%f,%f\n',range(temp+mm),aoa(temp+mm),vel(temp+mm));
               end
            end
        end
    end
   
    if (0)
        
    for ii=1:6  %every time write 6 groups of numdataPerObj data
        for nn=1:jj  % the tracker index
            temp = (nn-1)* numdataPerObj * 6 +1 + (ii-1)*numdataPerObj
            for mm=0:numdataPerObj-1
               % temp+mm
               if (genVarFlag==1)
                   fprintf(filePtr,'%f,%f,%f,%f,%f,%f\n',rangeT(temp+mm),aoaT(temp+mm),vel(temp+mm),varr,vara,varv);
               else
                   fprintf(filePtr,'%f,%f,%f\n',rangeT(temp+mm),aoaT(temp+mm),vel(temp+mm));
               end
            end
        end
    end
    end
end



%single tracker, multiple observations for allocation and association
%benchmark
index = zeros(1,4*200);
%index(201:400)=1;
%index(401:600)=2;
%index(601:end)=3;
symbolTemp = qammod(index,4) * 5;
symbol = awgn(symbolTemp,10);
xsymbol = real(symbol);
ysymbol = imag(symbol);

nextSymbolTemp = qammod(index,4,5/180*pi)*5;
nextsymbol = awgn(nextSymbolTemp,10);
nextxsymbol = real(nextsymbol);
nextysymbol = imag(nextsymbol);
range =[];
aoa = [];
vel = [];
fprintf(filePtr,'NUM=%d\n',800);
for jj=1:4*200
   range(jj) = sqrt(xsymbol(jj)^2+ysymbol(jj)^2);
   if (ysymbol(jj)>0)
      aoa(jj) = atan(xsymbol(jj)/ysymbol(jj));
   else
      aoa(jj) = atan(xsymbol(jj)/ysymbol(jj))+pi;
   end
   velx = nextxsymbol(jj)-xsymbol(jj);
   vely = nextysymbol(jj)-ysymbol(jj);
   vel(jj) = (xsymbol(jj)*velx+ysymbol(jj)*vely)/range(jj);
end
varr = [];
vara = [];
varv = [];
varr = 4*sum((range-mean(range)).^2)/(800-1);
vara = 4*sum((aoa-mean(aoa)).^2)/(800-1);
varv = 4*sum((vel-mean(vel)).^2)/(800-1);
for jj=1:4*200
    if (genVarFlag==1)
        fprintf(filePtr,'%f,%f,%f,%f,%f,%f\n',range(jj),aoa(jj),vel(jj),varr,vara,varv);
    else
        fprintf(filePtr,'%f,%f,%f\n',range(jj),aoa(jj),vel(jj));
    end
end

% next frame
fprintf(filePtr,'NUM=%d\n',800);
for jj=1:4*200
   range(jj) = sqrt(nextxsymbol(jj)^2+nextysymbol(jj)^2);
   if (nextysymbol(jj)>0)
      aoa(jj) = atan(nextxsymbol(jj)/nextysymbol(jj));
   else
      aoa(jj) = atan(nextxsymbol(jj)/nextysymbol(jj))+pi;
   end
   % speed no change, variance no change
   if (genVarFlag==1)
       fprintf(filePtr,'%f,%f,%f,%f,%f,%f\n',range(jj),aoa(jj),vel(jj),varr,vara,varv);
   else
       fprintf(filePtr,'%f,%f,%f\n',range(jj),aoa(jj),vel(jj));
   end

end
           
%figure;
%plot(real(symbol),imag(symbol),'.');
%hold on
%plot(real(nextsymbol),imag(nextsymbol),'ro');

% multiple trackers, for tracker update and prefict benchmark
index=[];
index = zeros(1,16*20);
for jj=1:15
    index(jj*20+1:(jj+1)*20)=jj;
end
symbolTemp = qammod(index,16) * 10;
symbol = awgn(symbolTemp,15);
xsymbol = real(symbol);
ysymbol = imag(symbol);

nextSymbolTemp = qammod(index,16,4/180*pi)*10;
nextsymbol = awgn(nextSymbolTemp,15);
nextxsymbol = real(nextsymbol);
nextysymbol = imag(nextsymbol);
%get first frame data
range =[];
aoa = [];
vel = [];
fprintf(filePtr,'NUM=%d\n',320);
for jj=1:320
   range(jj) = sqrt(xsymbol(jj)^2+ysymbol(jj)^2);
   if (ysymbol(jj)>0)
      aoa(jj) = atan(xsymbol(jj)/ysymbol(jj));
   else
      aoa(jj) = atan(xsymbol(jj)/ysymbol(jj))+pi;
   end
   velx = nextxsymbol(jj)-xsymbol(jj);
   vely = nextysymbol(jj)-ysymbol(jj);
   vel(jj) = (xsymbol(jj)*velx+ysymbol(jj)*vely)/range(jj);
end
varr = [];
vara = [];
varv = [];
varr = 4*sum((range-mean(range)).^2)/(320-1);
vara = 4*sum((aoa-mean(aoa)).^2)/(320-1);
varv = 4*sum((vel-mean(vel)).^2)/(320-1);
for jj=1:320
    if (genVarFlag==1)
        fprintf(filePtr,'%f,%f,%f,%f,%f,%f\n',range(jj),aoa(jj),vel(jj),varr,vara,varv); 
    else
        fprintf(filePtr,'%f,%f,%f\n',range(jj),aoa(jj),vel(jj));
    end
end

% next frame, increase the data number for update module benchmark
for ii=1:16  % 16 tracker
    numdataPerObj = ceil(100/ii);  % 100-600
    numdata = numdataPerObj *ii * 6;   % test 6 level of number of data
    index=[];
    index = zeros(1,numdata); 
    for jj=0:ii-1
         index(jj*(numdata/ii)+1:(jj+1)*(numdata/ii))=jj;
    end
    
    symbolTemp = qammod(index,16) * 10;
    symbolTempN = qammod(index,16,4/180*pi)*10;
    
    symbol = awgn(symbolTemp,15);
    xsymbol = real(symbol);
    ysymbol = imag(symbol);
    symbolN = awgn(symbolTempN,15);
    xsymbolN = real(symbolN);
    ysymbolN = imag(symbolN);
    range=[];
    aoa=[];
    vel=[];
    for jj=1:numdata
       range(jj) = sqrt(xsymbol(jj)^2+ysymbol(jj)^2);
       if (ysymbol(jj)>0)
          aoa(jj) = atan(xsymbol(jj)/ysymbol(jj));
       else
          aoa(jj) = atan(xsymbol(jj)/ysymbol(jj))+pi;
       end
       
       velx = xsymbolN(jj)-xsymbol(jj);
       vely = ysymbolN(jj)-ysymbol(jj);
       vel(jj) = (xsymbol(jj)*velx+ysymbol(jj)*vely)/range(jj);

    end
    varr = [];
    vara = [];
    varv = [];
    varr = 4*sum((range-mean(range)).^2)/(numdata-1);
    vara = 4*sum((aoa-mean(aoa)).^2)/(numdata-1);
    varv = 4*sum((vel-mean(vel)).^2)/(numdata-1);
    %numdata
    fprintf(filePtr,'NUM=%d\n',numdata);
    for jj=1:6  %every time write 6 groups of numdataPerObj data
        for nn=1:ii  % the tracker index
            temp = (nn-1)* numdataPerObj * 6 +1 + (jj-1)*numdataPerObj;
            for mm=0:numdataPerObj-1
               % temp+mm
               if (genVarFlag==1)
                   fprintf(filePtr,'%f,%f,%f,%f,%f,%f\n',range(temp+mm),aoa(temp+mm),vel(temp+mm),varr,vara,varv);
               else
                   fprintf(filePtr,'%f,%f,%f\n',range(temp+mm),aoa(temp+mm),vel(temp+mm));
               end
            end
        end
    end
end         






fclose(filePtr);



if 0
figure;
plot(squeeze(measureData(1,2,:))/3.1415926*180,'b.');
hold on
for ii=2:numFrame
    if (mod(ii,2)==0)
        plot(squeeze(measureData(ii,2,:))/3.1415926*180,'r.')
    else
        plot(squeeze(measureData(ii,2,:))/3.1415926*180,'b.')
    end
end


figure;
plot(squeeze(measureData(1,1,:)),'b.');
hold on
for ii=2:numFrame
    if (mod(ii,2)==0)
        plot(squeeze(measureData(ii,1,:)),'r.')
    else
        plot(squeeze(measureData(ii,1,:)),'b.')
    end
end

end
