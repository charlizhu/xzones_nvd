function rvOut = unrollRadialVelocity(rvMax, rvExp, rvIn)
% The distance between expected and measured values
if(rvExp >=0 )
    distance = rvExp - rvIn + rvMax;
    direction = sign(distance); 
else
    distance = rvExp - rvIn - rvMax;
    direction = sign(distance); 
end
factor = floor(abs(distance)/(2*rvMax));
rvOut = 2*rvMax*direction.*factor+rvIn;
end
