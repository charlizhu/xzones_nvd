classdef pointCloud < matlab.System
    %pointCloud
    %   Detailed explanation goes here
    
    properties (Nontunable)
        cloud
    end
    properties (Access = private, Nontunable)
        radarPos
        rangeMin
        rangeMax
        rangeStep
        numPoints
        xStd
        yCdf
        yVal
        % Velocity resolution of doppler FFT bin
        dopplerFftBinRes
        % Unambigious Maximum Radial Velocity, in +/- m/s
        uRadVelocityMax
    end
    % public methods
    methods
        % constructor
        function obj = pointCloud(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end
    
    methods (Access = protected)
        function resetImpl(obj)
            obj.radarPos = [0,0];
        end
        
        function setupImpl(obj)
            obj.rangeMin = obj.cloud.stat.ranges(1);
            obj.rangeMax = obj.cloud.stat.ranges(end);
            obj.rangeStep = obj.cloud.stat.ranges(2)-obj.cloud.stat.ranges(1);

            obj.numPoints.fRange = obj.cloud.stat.numPoints.fRange;
            obj.numPoints.velCoef = obj.cloud.stat.numPoints.velCoef;
            obj.numPoints.noise = 20*log10(obj.cloud.stat.numPoints.std);
            
            obj.dopplerFftBinRes = obj.cloud.dopplerFftBinRes;
            obj.uRadVelocityMax = obj.cloud.uRadVelocityMax;

            obj.xStd = obj.cloud.stat.xStd;
            obj.yCdf = obj.cloud.stat.ycdf;
            obj.yVal = obj.cloud.stat.yval;
        end
        
        function pCloud = stepImpl(obj, rec, vel)
            pCloud.xy = [];
            pCloud.u = [];
            pCloud.uvar = [];
            
            angleRes = 0.2; %in degrees
            angleResRad = angleRes*pi/180;
            rangeRes = obj.rangeStep; % in meters
            % Compute representation point
            sPoint = [rec(1)+rec(3)/2, rec(2)+rec(4)/2];
            % sRcsFactor = sqrt(rec(3)*rec(4)/(1.75*4)); %weight factor based on reflecting area ratio vs. standard car (1.75x4 m2)
            sRcsFactor = 1;
            % Compute range and angle
            sRange = sqrt((sPoint(1)-obj.cloud.radarPos(1))^2+(sPoint(2)-obj.cloud.radarPos(2))^2);
            sAngle = atan((sPoint(1)-obj.cloud.radarPos(1))/(sPoint(2)-obj.cloud.radarPos(2)));
            % Compute radial velocity
            sVel = vel(1)*sin(sAngle)+ vel(2)*cos(sAngle);
%{            
            if abs(sVel) < obj.dopplerFftBinRes
                return;
            end
%}
            % Compute number of points
            if(sRange < obj.rangeMin)
                rangeInd = 1;
            elseif (sRange > obj.rangeMax)
                rangeInd = length(obj.numPoints.fRange);
            else
                rangeInd = ceil((sRange-obj.rangeMin)/obj.rangeStep);
            end
            
            numRange = max(0, ceil(obj.numPoints.fRange(rangeInd)));
            numRange = ceil(numRange*sRcsFactor);
            if(numRange > 0)
                if(abs(sVel) < 3*obj.dopplerFftBinRes)
                    numSkirt = 0;
                else
                    % Compute "velocity skirt" factor. It is zero when radial
                    % velocity is exactly at half bin resolution, maximum +0.5
                    % or -0.5 at velocity bin edges
                    velSkirtFactor = (mod(sVel,obj.dopplerFftBinRes)-obj.dopplerFftBinRes/2)/obj.dopplerFftBinRes;
                    numSkirt = max(0, ceil(numRange*obj.numPoints.velCoef*abs(velSkirtFactor)));
                end
                numSkirtDist = [numSkirt-floor(numSkirt/3)-floor(numSkirt/4), floor(numSkirt/3), floor(numSkirt/4)];
                num = numRange + sum(numSkirtDist);
                
%                xstd = randn(1,1)*obj.xStd.std + obj.xStd.mean(rangeInd);
                xstd = obj.xStd.mean(rangeInd)/1.75*rec(3);
                x = randn(1,num)*xstd + sPoint(1);
                y = interp1(obj.yCdf, obj.yVal/4*rec(4), rand(1,num)) + sPoint(2);
                
                r = round(sqrt(x.^2+y.^2)/rangeRes)*rangeRes;
                a = round(atan(x./y)/angleResRad)*angleResRad;
                veld = round(sVel/obj.dopplerFftBinRes)*obj.dopplerFftBinRes;
                v = -veld*ones(1, numRange);
                for n=1:length(numSkirtDist)
                    if(numSkirtDist(n) > 0)
                        vel = veld+sign(velSkirtFactor)*n*obj.dopplerFftBinRes;
                        v = [v,-vel*ones(1,numSkirtDist(n))];
                    end
                end
                % Roll radial velocity into +/- uRadVelocityMax space
                rv=mod(v+obj.uRadVelocityMax,2*obj.uRadVelocityMax)-obj.uRadVelocityMax;
                
                %rvar = (r-sRange).^2;
                %avar = (a-sAngle).^2;
                %vvar = (v+sVel).^2;
                
                rvar = mean((r-sRange).^2)*ones(1, num);
                avar = mean((a-sAngle).^2)*ones(1, num);
                vvar = (obj.dopplerFftBinRes/2)^2*ones(1, num);
                
                x = r.*sin(a);
                y = r.*cos(a);
                
                pCloud.u = [r; a; rv];
                pCloud.uvar = [rvar; avar; vvar];
                pCloud.xy = [x; y];
            end
        end
    end    
end

