function relPoints = getCirclePart(nx, ny, xCenter, yCenter,rhoMin, rhoMax, minPhi, maxPhi)
% GETCIRCLEPART  Finds and returns the indices of points in a matrix
% corresponding to a desired part of a circle around a central position
%   nx, ny:  The x and y dimensions of the matrix for which the indices are
%   needed
%   xCenter, yCenter: The coordinates of the central position of the
%   particle
%   rhoMin, rhoMax: The minimum and maximum radii from the center for which
%   to include points
%   minPhi, maxPhi: The minimum and maximum angles (in radians) for which
%   to include points

[X,Y]=meshgrid(1:nx,1:ny);
c = [xCenter, yCenter];
angle = mod(atan2(Y-c(2),X-c(1)),2*pi);
minPhi = wrapTo2Pi(minPhi);
maxPhi = wrapTo2Pi(maxPhi);

% Checking which part of the circle to return (the smaller or larger part
% defined by the two angles)
if maxPhi > minPhi
    relPoints=  ((sqrt((X-c(1)).^2+(Y-c(2)).^2) >= rhoMin) & (sqrt((X-c(1)).^2+(Y-c(2)).^2) <= rhoMax)) & ...
                ((angle > minPhi) & (angle < maxPhi));
else
    relPoints=  ((sqrt((X-c(1)).^2+(Y-c(2)).^2) >= rhoMin) & (sqrt((X-c(1)).^2+(Y-c(2)).^2) <= rhoMax)) & ...
                ((angle > minPhi) | (angle < maxPhi));    
end