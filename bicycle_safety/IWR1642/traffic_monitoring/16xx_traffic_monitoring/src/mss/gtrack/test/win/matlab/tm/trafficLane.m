classdef trafficLane < matlab.System
    properties (Nontunable)
        %lane parameters
        lane
        %handle to pCloud object
        pCloudHandle
        % Axes
        ax
    end
    properties (Access = private, Nontunable)
        pMaxNumVehicles = 100 %Maximum number of active vehicles per lane
        pMaxNumVehiclesIDs = 1024
        pXpos
        pWidth
        pEnterLine
        pStopLine
        pExitLine
        pNumVehicles
        pVehicleHandle
        pVehicleFreeList
        pVehicleActiveList
        pHeartBeatCount
        pRecHandle
        pPlotHandle
        pColorArray
    end
    
    % public methods
    methods
        % constructor
        function obj = trafficLane(varargin)
            setProperties(obj, nargin, varargin{:});
        end
        function hVeh = getVehicleHandle(obj, id)
            hVeh = obj.pVehicleHandle{id};
        end
       function setVehicleMark(obj, id, marker)
           obj.pVehicleHandle{id}.marker = marker;
       end
       function marker = getVehicleMark(obj, id)
           if isa(obj.pVehicleHandle{id},'handle') && isvalid(obj.pVehicleHandle{id})
               marker = obj.pVehicleHandle{id}.marker;
           else
               marker = -1;
           end
       end
    end
    
    methods (Access = protected)
        function resetImpl(obj)
            obj.pHeartBeatCount = 0;
        end
        function setupImpl(obj)
            obj.pWidth = obj.lane.width;
            obj.pXpos = obj.lane.xpos;
            obj.pEnterLine = obj.lane.enterline;
            obj.pStopLine = obj.lane.stopline;
            obj.pExitLine = obj.lane.exitline;
            
            obj.pVehicleHandle = cell(obj.pMaxNumVehicles,1);
            obj.pVehicleFreeList = 1:obj.pMaxNumVehiclesIDs;
            obj.pVehicleActiveList =[];
            obj.pRecHandle = zeros(obj.pMaxNumVehicles,1);
            obj.pPlotHandle = zeros(obj.pMaxNumVehicles,1);
            obj.pColorArray = 'ymcrgbk';
        end
        
        function laneStat = stepImpl(obj, arrivals, light, dt)
            obj.pHeartBeatCount  = obj.pHeartBeatCount +1;
            if(~strcmp(arrivals.type,'none'))
                % we have new vehicle arrival

                vid = allocateVehicle(obj);                 % allocate new Vehicle ID
                hid = mod(vid-1,obj.pMaxNumVehicles) +1;    % we have only 100 handles
                pos = [obj.pXpos+obj.pWidth/2 obj.pEnterLine]; % [posX posY]
                vel = [0 arrivals.speed]; %[velX velY];
                obj.pVehicleHandle{hid} = vehicle('entryTime', arrivals.tick, 'id', vid, 'lane', obj.lane.num, 'type', arrivals.type, 'velocity', vel, 'position', pos);
                disp([num2str(obj.pHeartBeatCount),': ENTER vid=',num2str(obj.lane.num*obj.pMaxNumVehiclesIDs+vid),', type=',arrivals.type,', pos=(',num2str(pos(1),2),',',num2str(pos(2),2),'), vel=',num2str(abs(vel(2)),2)]);
            end
            
            activeList = obj.pVehicleActiveList;
            movingVehicleCount = 0;
            pendingVehicleCount = 0;

            pCloudLane.id = [];
            pCloudLane.xy = [];
            pCloudLane.u = [];
            pCloudLane.uvar = [];
            
            % determine the first vehicle behaviour
            if(~isempty(activeList))
                vid = activeList(1);
                hid = mod(vid-1,obj.pMaxNumVehicles) + 1;
                posy = obj.pVehicleHandle{hid}.position(2);
                vely = obj.pVehicleHandle{hid}.velocity(2);

                if(posy < obj.pStopLine)
                    light = 'Green';
                end
                switch light
                    case 'Yellow'
                        spaceNeeded = (vely^2)/(2*10);
                        if(posy > spaceNeeded)
                            obstacle.Line = obj.pStopLine;
                            obstacle.Speed = 0;
                        else
                            obstacle.Line = -100;
                            obstacle.Speed = 0;
                        end
                    case 'Red'
                        obstacle.Line = obj.pStopLine;
                        obstacle.Speed = 0;
                    otherwise
                        obstacle.Line = -100;
                        obstacle.Speed = 0;
                end
 
                for n = 1:length(activeList)
                    vid = activeList(n);
                    hid = mod(vid-1,obj.pMaxNumVehicles)+1;
                    % Move the vehicle
                    [rec, vel] = step(obj.pVehicleHandle{hid}, obstacle, dt);
                    obstacle.Line = rec(2) + rec(4) + 1;
                    obstacle.Speed = vel(2);

                    % Generate Point Cloud
                    pCloudVeh = step(obj.pCloudHandle, rec, vel);
                    
                    if(obj.pPlotHandle(hid)~=0)
                        delete(obj.pPlotHandle(hid));
                        obj.pPlotHandle(hid) = 0;
                    end
                    if ~isempty(pCloudVeh.xy)
                        obj.pPlotHandle(hid) = plot(obj.ax(2), pCloudVeh.xy(1,:), pCloudVeh.xy(2,:),'.', 'Color', obj.pColorArray(mod(hid,length(obj.pColorArray))+1));
                    end
                    
                    if(rec(2) < 10)
                        disp([num2str(obj.pHeartBeatCount),': EXIT vid=',num2str(obj.lane.num*obj.pMaxNumVehiclesIDs+vid)]);
                        setVehicleMark(obj, hid, 1); % Mark the vehicle data as ready for logging
                        freeVehicle(obj, vid);
                        obj.pVehicleHandle{hid}.state = 0;
%                        delete(obj.pVehicleHandle{id});
                        delete(obj.pRecHandle(hid));
                        delete(obj.pPlotHandle(hid));
                        obj.pRecHandle(hid) = 0;
                        obj.pPlotHandle(hid) = 0;
                    else
                        if(obj.pRecHandle(hid)~=0)
                            delete(obj.pRecHandle(hid));
                        end
                        if(obj.pVehicleHandle{hid}.state == 1)
                            pendingVehicleCount = pendingVehicleCount +1;
                        else
                            obj.pRecHandle(hid) = rectangle(obj.ax(1), 'Position', rec, 'EdgeColor', 'k');
                            if vel(2) > 0.5
                                movingVehicleCount = movingVehicleCount +1;
                            end
                        end
                    end
                    if ~isempty(pCloudVeh.xy)
                        idArray = ones(1,size(pCloudVeh.xy,2))*obj.lane.num*obj.pMaxNumVehiclesIDs+vid;
                        pCloudLane.id = [pCloudLane.id idArray];
                        pCloudLane.xy = [pCloudLane.xy pCloudVeh.xy];
                        pCloudLane.u = [pCloudLane.u pCloudVeh.u];
                        pCloudLane.uvar = [pCloudLane.uvar pCloudVeh.uvar];
                    end
                end
            end
            laneStat.total = length(obj.pVehicleActiveList);
            laneStat.moving = movingVehicleCount;
            laneStat.pending = pendingVehicleCount;
            laneStat.pcloud = pCloudLane;
        end
        
        function id = allocateVehicle(obj)
            %allocateVehicle Allocate new vehicle ID from the free pool
            % Get new vehicle ID
            [id, obj.pVehicleFreeList] = getFromList(obj.pVehicleFreeList);
            obj.pNumVehicles = obj.pNumVehicles + 1;
            % Adds vehicle ID to the active list
            obj.pVehicleActiveList = addToList(obj.pVehicleActiveList, id);
        end
        
        function freeVehicle(obj, id)
            %freeVehicle Places vehicle ID back to free pool
            % Remove vehicle ID from the active list
            obj.pVehicleActiveList = removeFromList(obj.pVehicleActiveList, id);
            % Returns vehicle id back to the free pool
            obj.pVehicleFreeList = addToList(obj.pVehicleFreeList, id);
            obj.pNumVehicles = obj.pNumVehicles - 1;
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