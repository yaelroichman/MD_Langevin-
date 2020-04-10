function [fx,fy]=getWCAForces(particlesMat, R, epsilon, particleRepel, wallRepel, wallPositionsX, wallPositionsY)
% particlesMat: The positions of the particles in the current step
% R: The radii of the different particles (Only works for a single R
% currently)
% wallRepel: A boolean to determine whether or not to compute repulsion from
% walls
% particleRepel:  A boolean to determine whether or not to compute
% repulsion between the particles
% wallPositionX: The x coordinates of the walls
% wallPositionY: The y coordinates of the walls
% In WCA notation I take sigma->R
numOfParticles = size(particlesMat,1);
fx = zeros(numOfParticles,1);
fy = zeros(numOfParticles,1);
Rc = 2.5*R;

for currParticle = 1:numOfParticles
    checkedParticlePosition = particlesMat(currParticle,:);
    if wallRepel
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
    %% Checking interactions with other particles
    if particleRepel
        for pairedParticle = (currParticle+1):numOfParticles
            pairedParticlePosition = particlesMat(pairedParticle,:);
            yDiff = checkedParticlePosition(2)-pairedParticlePosition(2);
            xDiff = checkedParticlePosition(1)-pairedParticlePosition(1);

            distance = sqrt(xDiff^2+yDiff^2) - R;
            if distance <= Rc
                angle = atan2(yDiff,...
                              xDiff);
                totalForce = getLJForce(distance, epsilon, R);

                fx(currParticle) = fx(currParticle) + totalForce.*cos(angle);
                fx(pairedParticle) = fx(pairedParticle) - totalForce.*cos(angle);
                fy(currParticle) = fy(currParticle) + totalForce.*sin(angle);
                fy(pairedParticle) = fy(pairedParticle) - totalForce.*sin(angle);
            end
        end
    end
end
