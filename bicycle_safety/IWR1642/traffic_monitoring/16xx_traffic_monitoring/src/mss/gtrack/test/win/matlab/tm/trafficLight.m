classdef trafficLight < matlab.System
    properties (Nontunable)
        %lane parameters
    timeouts
    end
properties (Constant = true)
    STATE_GREEN = 1;
    STATE_YELLOW = 2;
    STATE_RED = 3;
end     
    properties (Access = private, Nontunable)
        pState
        pActiveTime
        pRedTime
        pYellowTime
        pGreenTime
    end
    % public methods
    methods
        % constructor
        function obj = trafficLight(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end
    
    methods (Access = protected)
        function resetImpl(obj)
            obj.pState = obj.STATE_GREEN;
            obj.pActiveTime = 0;
        end
        function setupImpl(obj)
            obj.pGreenTime = obj.timeouts(1);
            obj.pYellowTime = obj.timeouts(2);
            obj.pRedTime = obj.timeouts(3);
        end
        
        function light = stepImpl(obj, dt)
            obj.pActiveTime = obj.pActiveTime + dt;
            switch obj.pState
                case obj.STATE_GREEN
                    if obj.pActiveTime >= obj.pGreenTime
                        if(obj.pYellowTime ~=0)
                            obj.pState = obj.STATE_YELLOW;
                        elseif(obj.pRedTime ~=0)
                            obj.pState = obj.STATE_RED;
                        end
                        obj.pActiveTime = 0;
                    end
                case obj.STATE_YELLOW
                    if obj.pActiveTime >= obj.pYellowTime
                        obj.pState = obj.STATE_RED;
                        obj.pActiveTime = 0;
                    end
                case obj.STATE_RED
                    if obj.pActiveTime >= obj.pRedTime
                        obj.pState = obj.STATE_GREEN;
                        obj.pActiveTime = 0;
                    end
            end
            light = obj.pState;
        end
    end
end