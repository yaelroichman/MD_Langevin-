function total_force = getTrapForcesForSingleParticle(particlePosition, trapPositions, A, s)
%Calculate force from all optical traps on a particle
f = (A./s.^2) .* (particlePosition-trapPositions) .* exp((-0.5*(particlePosition-trapPositions)./s).^2);
total_force = sum(f,1);