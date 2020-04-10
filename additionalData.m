classdef additionalData < handle
    properties
        wallMoveSteps % The step numbers in which the wall was moved
        closestParticlePositions % The positions of the closest particle in each wall move
        newWallPositions % The positions to which the wall was moved
        wallMoveInd;
    end
end