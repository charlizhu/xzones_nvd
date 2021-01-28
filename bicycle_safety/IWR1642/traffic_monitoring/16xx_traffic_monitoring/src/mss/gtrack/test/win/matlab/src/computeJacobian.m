function [J] = computeJacobian(format, s)

% Function to compute the Jacobian
% Muhammad Ikram - mzi@ti.com; June, 2013
% Copyright - Texas Instruments Incorporated

% 0: (2D) Input s = [posx posy velx vely],
%       Output: J = 3x4 
% 1: (2DA) Input s = [posx posy velx vely 0 0],
%       Output: J = 3x6 
% 2: (3D) Input s = [posx posy posz velx vely velx] 
%       Output: J = 4x6
if(format == 0)
    J = zeros(3,4);
else
    J = zeros(3,6);
end    
switch format
    case{0,1}
        posx = s(1); posy = s(2); velx = s(3); vely = s(4);
        r = sqrt(posx^2 + posy^2);
        J(1,1) = posx/r;
        J(1,2) = posy/r;
        
        J(2,1) = posy/(posx^2 + posy^2);
        J(2,2) = -posx/(posx^2 + posy^2);
        
        J(3,1) = (posy*(velx*posy - posx*vely))/r^3;
        J(3,2) = (posx*(posx*vely - velx*posy))/r^3;
        J(3,3) = posx/r;
        J(3,4) = posy/r;
    otherwise
        % TODO
        error('radarTracking:invalidParameter', 'CoputeJformatNotSupported');
end
end