function total_force = getTrapForcesForSingleParticle(particlePosition, trapPositions, A, s)
%Calculate force from all optical traps on a particle
%f = (A./s.^2) .* (particlePosition-trapPositions) .* exp(-((0.5*(particlePosition-trapPositions)./s).^2));
%total_force = sum(f,1);
%%
f_x=(A(:,1)./s(:,1).^2).*(particlePosition(:,1)-trapPositions(:,1)).*exp(-((0.5*(particlePosition(:,1)-trapPositions(:,1))./s(:,1)).^2+(0.5*(particlePosition(:,2)-trapPositions(:,2))./s(:,2)).^2));
f_y=(A(:,2)./s(:,2).^2).*(particlePosition(:,2)-trapPositions(:,2)).*exp(-((0.5*(particlePosition(:,1)-trapPositions(:,1))./s(:,1)).^2+(0.5*(particlePosition(:,2)-trapPositions(:,2))./s(:,2)).^2));
f=[f_x,f_y];
total_force = sum(f,1);
