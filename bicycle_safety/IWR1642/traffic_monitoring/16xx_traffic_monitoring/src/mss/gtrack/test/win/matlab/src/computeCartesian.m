function [cart] = computeCartesian(format, sph)

% Function to compute the non-linear conversion from Spherical to Cartesian
% coordinates 
% Copyright - Texas Instruments Incorporated 
% Input formats supported:
% Input formats supported:
% 0: 2D s = [posx posy velx vely] or
% 1: 2DA s = [posx posy velx vely accx accy]
% 2: 3D s = [posx posy posz velx vely velx]

% 0: (2D) Input sph = [range azimuth doppler],
%       Output: cart = [posx posy velx vely]T 
% 1: (2DA) Input sph = [range azimuth doppler],
%       Output: cart = [posx posy velx vely 0 0]T 
% 2: (3D) Input sph = [range elevation azimuth doppler] 
%       Output: cart = [posx posy posz velx vely velx]T

switch format
    case {0,1}
        r = sph(1); phi = sph(2); vr = sph(3);
        posx = r*sin(phi);
        posy = r*cos(phi);

        velx = vr*sin(phi);
        vely = vr*cos(phi);

        if(format == 0)
            cart = [posx, posy, velx, vely]';
        else
            cart = [posx, posy, velx, vely, 0, 0]';
        end
        
    case 2
        r = sph(1); phi = sph(2); omega = sph(3); vr = sph(4);
        posx = r*cos(phi)*sin(omega);
        posy = r*sin(phi)*sin(omega);
        posz = r*cos(omega);
        
        velx = vr*cos(phi)*sin(omega);
        vely = vr*sin(phi)*sin(omega);
        velz = vr*cos(omega);
        cart = [posx, posy, posz, velx, vely, velz]';
end
end