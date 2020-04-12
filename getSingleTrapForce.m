function f = getSingleTrapForce(particlePosition, trapPosition, A, s)
%Calculate force from an optical trap on a particle
f = (A/s^2) * (particlePosition-trapPosition) * exp(-0.5*(particlePosition-trapPosition)/s)^2;