function [particlesX,particlesY] = randomizePositions(wallPositionsX, wallPositionsY, numOfParticles, R)
% Randomizes the particle positions while ensuring no particles overlap
particlesX = zeros(numOfParticles,1);
particlesY = zeros(numOfParticles,1);
rng(0,'twister');
lx = abs(wallPositionsX(2) - wallPositionsX(1));
ly = abs(wallPositionsY(2) - wallPositionsY(1));
for currParticle = 1:numOfParticles
    collision = true;
    %Continuing until there is no collision
    while collision
        collision = false;
        % Randomizing a position for the current particle
        particlesX(currParticle) = (lx-2.*R).*rand(1,1) + wallPositionsX(1) + R;
        particlesY(currParticle) = (ly-2.*R).*rand(1,1) + wallPositionsY(1) + R;
        for collidedParticle = 1:currParticle
            currCollision = ...
                checkCollision(particlesX(currParticle), particlesY(currParticle), R,...
                particlesX(collidedParticle),particlesY(collidedParticle), R);
            if currCollision
                collision = true;
                break
            end
        end
    end
end