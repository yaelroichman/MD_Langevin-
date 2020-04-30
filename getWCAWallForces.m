function [fx,fy]=getWCAWallForces(particlesMat, R, epsilon, wallPositionsX, wallPositionsY)
% getWCAWallForces computes the forces on the particles from the walls
% using WCA repulsion
% particlesMat: The positions of the particles in the current step
% R: The radii of the different particles (Only works for a single R
% currently)
% wallPositionX: The x coordinates of the walls
% wallPositionY: The y coordinates of the walls
numOfParticles = size(particlesMat,1);
fx = zeros(numOfParticles,1);
fy = zeros(numOfParticles,1);
Rc = 2.5*R;

for currParticle = 1:numOfParticles
    checkedParticlePosition = particlesMat(currParticle,:);
    distanceFromWallsX = abs(checkedParticlePosition(1).*ones(1,2) - wallPositionsX);
    distanceFromWallsY = abs(checkedParticlePosition(2).*ones(1,2) - wallPositionsY);
    
    %% Checking if the particle is near a wall
    if distanceFromWallsX(1) <= Rc
        fx(currParticle) = fx(currParticle) + getLJForce(distanceFromWallsX(1),epsilon, R);    
    end
    if distanceFromWallsX(2) <= Rc
        fx(currParticle) = fx(currParticle) - getLJForce(distanceFromWallsX(2),epsilon, R);
    end
    if distanceFromWallsY(1) <= Rc
        fy(currParticle) = fy(currParticle) + getLJForce(distanceFromWallsY(1),epsilon, R);    
    end
    if distanceFromWallsY(2) <= Rc
        fy(currParticle) = fy(currParticle) - getLJForce(distanceFromWallsY(2),epsilon, R);
    end
end