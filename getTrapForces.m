function [fx,fy]=getTrapForces(particlePositions, trapPositions, A, s)
%%Calculate the forces all of the traps apply on all of the particles.

numOfParticles = size(particlePositions,1);
numOfTraps = size(trapPositions,1);

%If A,s (depth, size of the trap) are scalars convert them into a matrix.
if numel(A) == 1
    A = A.*ones(numOfTraps,2);
end
if numel(s) == 1
    s = s.*ones(numOfTraps,2);
end

fx = zeros(numOfParticles,1);
fy = zeros(numOfParticles,1);

for currParticle = 1:numOfParticles
    checkedParticlePosition = particlePositions(currParticle,:);
    for currTrap = 1:numOfTraps
        checkedTrapPosition = trapPositions(currTrap,:);
        checkedTrapDepth = A(currTrap,:);
        checkedTrapSize = s(currTrap,:);
        
        fx(currParticle) = fx(currParticle) + getSingleTrapForce(checkedParticlePosition(1),...
            checkedTrapPosition(1), checkedTrapDepth(1), checkedTrapSize(1));
        
        fy(currParticle) = fy(currParticle) + getSingleTrapForce(checkedParticlePosition(2),...
            checkedTrapPosition(2), checkedTrapDepth(2), checkedTrapSize(2));
        
    end
end