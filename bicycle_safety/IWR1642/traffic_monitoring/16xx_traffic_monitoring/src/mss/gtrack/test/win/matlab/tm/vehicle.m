classdef vehicle < matlab.System
   properties
       %entry time stamp
       entryTime
       % vehicle ID 
       id
       % Index of the lane the vehicle in {1:6} 
       lane
       % Vehicle type {truck | car | motorcycle | bicylce}
       type = 'car'
       % Vehicle states = {stopped at traffic light | constant speed approaching | accelerating | decelerating | manuvering right | manuevering left} 
       state
       % Vehicle geometry widht, height, length
       width
       height
       length
       % Vehicle position, relative to the stop line
       position = [0,0]
       % vehicle velocity: (velX (horizontal, lane changing), velY (vertical))
       velocity = [0,0]
       % vehicle acceleration: (accX, accY)
       acceleration = [0,0]
       %rectangle
       rectangle = [0,0,0,0]
       %maxVelocity
       maxVelocity
       %maxAcc
       maxAccel
       %marker
       marker
       % vehicle history
       tick
       posxyHistory
       velxyHistory
       accxyHistory
   end
properties (Constant = true)
    STATE_PENDING = 1;
    STATE_VCONST = 2;
    STATE_DECEL = 3;
    STATE_STOPPED = 4;
    STATE_ACCEL = 5;
    STATE_LEFT = 6;
end   
   % public methods
   methods
       % constructor
       function obj = vehicle(varargin)
           setProperties(obj, nargin, varargin{:});
           obj.state = obj.STATE_PENDING;
           obj.acceleration = [0 0];
           obj.posxyHistory = zeros(1000,2);
           obj.velxyHistory = zeros(1000,2);
           obj.accxyHistory = zeros(1000,2);
           obj.tick = 1;
           obj.marker = 0;
           switch obj.type
               case 'car'
                   obj.width = 1.75;
                   obj.height = 1.7;
                   obj.length = 4;
                   obj.maxVelocity = 40; %max velocity, m/s (40m/s = 144km/h)
                   obj.maxAccel = 10; %max acceleration, m/s2 (ex. with tires static friction of 0.75 => max acc = 7.35m/s2)
               case 'truck'
                   obj.width = 2.55;
                   obj.height = 1.4;
                   obj.length = 8;
                   obj.maxVelocity = 25; %max velocity, m/s (25m/s = 90km/h)
                   obj.maxAccel = 4; %max acceleration, m/s2 (2m/s2)
               case 'motorcycle'
                   obj.width = 1.0;
                   obj.height = 1.4;
                   obj.length = 3;
                   obj.maxVelocity = 40; %max velocity, m/s (40m/s = 144km/h)
                   obj.maxAccel = 10; %max acceleration, m/s2 (4m/s2)
                   
               case 'bicycle'
                   obj.width = 0.5;
                   obj.height = 1.4;
                   obj.length = 2;
                   obj.maxVelocity = 8; %max velocity, m/s (8m/s = 29km/h)
                   obj.maxAccel = 5; %max acceleration, m/s2 (4m/s2)
           end
           
           obj.rectangle = [-obj.width/2 0 obj.width obj.length];
           obj.posxyHistory(1,:) = obj.position;
           obj.velxyHistory(1,:) = obj.velocity;
       end
   end
   methods (Access = protected)
       function resetImpl(obj)
       end
       
       function setupImpl(obj)
       end
       
       function [rec, vel] = stepImpl(obj, obstacle, dt)           
           switch obj.state
               case obj.STATE_PENDING
                   if(obj.position(2) > obstacle.Line + 2)
                       if((obstacle.Line > 0) && (obj.velocity(2) > obstacle.Speed))
                           %if there is an obstacle in front, compute the space needed
                           spaceNeeded = (obj.velocity(2)^2-obstacle.Speed^2)/(2*obj.maxAccel);
                           if((obj.position(2) - obstacle.Line) > spaceNeeded + 10)
                               % wait till we have room to approach
                               obj.state = obj.STATE_VCONST;
                           end
                       else
                            obj.state = obj.STATE_VCONST;
                       end
                   end
                   if obj.state == obj.STATE_PENDING
                       obj.entryTime = obj.entryTime + 1;
                   end
                   
               case obj.STATE_VCONST
                   if((obstacle.Line > 0) && (obj.velocity(2) > obstacle.Speed))
                       % if we appraching an obstacle
                       spaceNeeded = (obj.velocity(2)^2-obstacle.Speed^2)/(2*obj.maxAccel);                       
                       if((obj.position(2) - obstacle.Line) < spaceNeeded)
                           obj.acceleration(2) = -obj.maxAccel;               
                           obj.state = obj.STATE_DECEL;
                       end
                   end
               case obj.STATE_DECEL
                   if((obstacle.Line > 0) && (obj.velocity(2) > obstacle.Speed))
                       % We are appraching an obstacle which is close enough, then update deceleration value
                       timeToDecel = 2*(obj.position(2) - obstacle.Line)/(obj.velocity(2) + obstacle.Speed);
                       obj.acceleration(2) = -(obj.velocity(2)-obstacle.Speed)/timeToDecel;
                   else
                       % no more obstacles, or the obstacle is faster, accelerate
                       obj.acceleration(2) = 0;
                       obj.state = obj.STATE_ACCEL;
                   end
                   
               case obj.STATE_STOPPED
                   if(obj.position(2) - obstacle.Line > 2)
                       obj.acceleration(2) = obj.maxAccel/8;
                       obj.state = obj.STATE_ACCEL;
                   end
                   
               case obj.STATE_ACCEL
                   if((obstacle.Line > 0) && (obj.velocity(2) > obstacle.Speed))
                       % if we approaching an obstacle
                       obj.state = obj.STATE_DECEL;
                       obj.acceleration(2) = -obj.maxAccel/8;
                   else
                       obj.acceleration(2) = min(obj.acceleration(2) + obj.maxAccel/8, obj.maxAccel);
                   end
           end
           
               
           if(obj.state ~= obj.STATE_PENDING)
               % update velocity
               newVelocity = obj.velocity + obj.acceleration*dt;
              
               % Don't allow to accelerate beyond max velocity
               if(obj.state == obj.STATE_ACCEL)
                   if(newVelocity(2) > obj.maxVelocity)
                       obj.acceleration(2) = 0;
                       obj.state = obj.STATE_VCONST;
                   end
               end
             
               if(newVelocity(2) <= 0.01)
                   obj.state = obj.STATE_STOPPED;
                   newVelocity(2) = 0;
                   obj.acceleration(2) = 0;
               end
               obj.velocity = newVelocity;
               
               
               % update position
               obj.position = obj.position - 0.5*(obj.velocity+newVelocity)*dt;
           end
           
           rec = [obj.position 0 0] + obj.rectangle;
           vel = obj.velocity;
           
           if(obj.state ~= obj.STATE_PENDING)
               ind = mod(obj.tick,1000)+1;
               obj.posxyHistory(ind,:) = obj.position;
               obj.velxyHistory(ind,:) = obj.velocity;
               obj.accxyHistory(ind,:) = obj.acceleration;
               obj.tick = obj.tick + 1;
           end
       end
   end
end