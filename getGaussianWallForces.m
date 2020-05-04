function [fx,fy]=getGaussianWallForces(particlesMat, R, A,s, wallPositionsX, wallPositionsY)
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
    xBeyondLeftWall = wallPositionsX(1) - (checkedParticlePosition(1) - R);
    xBeyondRightWall = checkedParticlePosition(1) + R - wallPositionsX(2);
    yBeyondTopWall = checkedParticlePosition(2) + R - wallPositionsY(2);
    yBeyondBottomWall = wallPositionsY(1) - (checkedParticlePosition(2) - R);
    
    %% Checking if the particle's edge has passed a wall
    if xBeyondRightWall > 0
        fx(currParticle) = fx(currParticle) - getGaussianForce(A,s,xBeyondRightWall);
    end
    if xBeyondLeftWall > 0
        fx(currParticle) = fx(currParticle) + getGaussianForce(A,s,xBeyondLeftWall);
    end
    if yBeyondTopWall > 0
        fy(currParticle) = fy(currParticle) - getGaussianForce(A,s,yBeyondTopWall);
    end
    if yBeyondBottomWall > 0
        fy(currParticle) = fy(currParticle) + getGaussianForce(A,s,yBeyondBottomWall);
    end
end