function [fx,fy]=getWCAParticleForces(particlesMat, R, epsilon)
% getWCAParticleForces computes the forces on the particles from the other
% particles using WCA repulsion (effectively similar to hard sphere)
% particlesMat: The positions of the particles in the current step
% R: The radii of the different particles (Only works for a single R
% currently)
numOfParticles = size(particlesMat,1);
fx = zeros(numOfParticles,1);
fy = zeros(numOfParticles,1);
Rc = 2.5*R;

for currParticle = 1:numOfParticles
    checkedParticlePosition = particlesMat(currParticle,:);
    
    %% For each particle, checking its interaction with all others
    for pairedParticle = (currParticle+1):numOfParticles
        pairedParticlePosition = particlesMat(pairedParticle,:);
        yDiff = checkedParticlePosition(2)-pairedParticlePosition(2);
        xDiff = checkedParticlePosition(1)-pairedParticlePosition(1);

        distance = sqrt(xDiff^2+yDiff^2) - R;
        if distance <= Rc
            %% Computing the forces between the particles
            totalForce = getLJForce(distance, epsilon, R);
            angle = atan2(yDiff,...
                          xDiff);

            %% Adding the force to both particles with opposite signs
            fx(currParticle) = fx(currParticle) + totalForce.*cos(angle);
            fx(pairedParticle) = fx(pairedParticle) - totalForce.*cos(angle);
            fy(currParticle) = fy(currParticle) + totalForce.*sin(angle);
            fy(pairedParticle) = fy(pairedParticle) - totalForce.*sin(angle);
        end
    end
end