function [H] = computeH(format, s)

% Function to compute the non-linear conversion from Cartesian to Spherical
% coordinates
% Muhammad Ikram - mzi@ti.com; June, 2013
% Copyright - Texas Instruments Incorporated

% Input formats supported:
% 0: 2D s = [posx posy velx vely] or
% 1: 2DA s = [posx posy velx vely accx accy]
% 2: 3D s = [posx posy posz velx vely velx]
% Output formats supported:
% H = [range azimuth doppler]T (2D or 2DA) or
% H = [range elevation azimuth doppler]T (3D)
switch format
    case {0,1}
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
    otherwise
        % TODO
        error('radarTracking:invalidParameter', 'ComputeHformatNotSupported');
end
return

