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
    xRelToLeftWall = abs(wallPositionsX(1) - checkedParticlePosition(1));
    xRelToRightWall = abs(checkedParticlePosition(1) - wallPositionsX(2));
    yRelToTopWall = abs(checkedParticlePosition(2) - wallPositionsY(2));
    yRelToBottomWall = abs(wallPositionsY(1) - checkedParticlePosition(2));
    

    fx(currParticle) = fx(currParticle) - getGaussianForce(A,s,xRelToRightWall);

    fx(currParticle) = fx(currParticle) + getGaussianForce(A,s,xRelToLeftWall);

    fy(currParticle) = fy(currParticle) - getGaussianForce(A,s,yRelToTopWall);

    fy(currParticle) = fy(currParticle) + getGaussianForce(A,s,yRelToBottomWall);

end
