function [XX, YY, ZZ, v] = gatePlot3(ax, G, C, A)

%Extract the ellipsoid's axes lengths (a,b,c) and the rotation matrix (V) using singular value decomposition:
[U,D,V] = svd(A/G);

a = 1/sqrt(D(1,1));
b = 1/sqrt(D(2,2));
c = 1/sqrt(D(3,3));
v = 4*pi*a*b*c/3;
    
% generate ellipsoid at 0 origin
[X,Y,Z] = ellipsoid(0,0,0,a,b,c);
XX = zeros(size(X));
YY = zeros(size(X));
ZZ = zeros(size(X));
for k = 1:length(X)
    for j = 1:length(X)
        point = [X(k,j) Y(k,j) Z(k,j)]';
        P = V * point;
        XX(k,j) = P(1)+C(1);
        YY(k,j) = P(2)+C(2);
        ZZ(k,j) = P(3)+C(3);
    end
end
%disp(max(max(YY)));
%{
x = XX.*sin(YY);
y = XX.*cos(YY);
mesh(ax, x,y,ZZ);
view(ax, [0,90]);
hidden off
%}