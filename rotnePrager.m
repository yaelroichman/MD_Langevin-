function DMat = rotnePrager(particlesMat,R, D)
d = 2; % The dimension of the problem. 
% The code is NOT scalable to more or less dimensions, the variable d is
% only present for better readability.

% Converting the matrix of particle positions into a vector of the form
% (x1, y1, x2, y2,..., xn, yn)
numOfParticles = size(particlesMat,1);
workingMat = particlesMat';
particlesVec = zeros(d*numOfParticles,1);
for i = 1:numOfParticles
    particlesVec(2*i-1) = workingMat(1,i);
    particlesVec(2*i) = workingMat(2,i);
end

% Constants for the matrix values
c1 = 0.75*R;
c2 = 0.5*R^3;

% Setting a unit matrix to start from
DMat = eye(d*numOfParticles);

for checkedParticle = 1:numOfParticles
    for pairedParticle = 1:(checkedParticle-1)
        % Getting the relevant matrix positions
        xxPosition = [checkedParticle*d-1, pairedParticle*d-1];
        xyPosition = [checkedParticle*d, pairedParticle*d-1];
        yxPosition = [checkedParticle*d-1, pairedParticle*d];
        yyPosition = [checkedParticle*d, pairedParticle*d];
        
        % Getting the distances in the axes, and in total
        xDiff = particlesMat(checkedParticle,1) - particlesMat(pairedParticle,1);
        yDiff = particlesMat(checkedParticle,2) - particlesMat(pairedParticle,2);
        r = sqrt(xDiff^2+yDiff^2);
        
        % For the x-x matrix position, adding c1(1 + x^2/r^2)
        DMat(xxPosition(1),xxPosition(2)) = (c1./r).*(1 + xDiff^2/r^2) +...
                                            (c2 ./ r^3).*(1-(xDiff^2/r^2));
        
        % For the y-y matrix position, adding c1(1 + y^2/r^2)
        DMat(yyPosition(1), yyPosition(2)) = (c1./r).*(1 + yDiff^2/r^2) +...
                                             (c2 ./ r^3).*(1-yDiff^2/r^2);
        
        % For the x-y position, adding c2(1-3xy/r^2)
        DMat(xyPosition(1), xyPosition(2)) = (c1./r).*((xDiff.*yDiff)./r^2) -...
                                             3.*(c2./r^3).*((xDiff.*yDiff./r.^2));
        DMat(yxPosition(1), yxPosition(2)) = DMat(xyPosition(1), xyPosition(2));
    end
end

DMat = DMat.*D;
DMat = DMat'+DMat-diag(diag(DMat));