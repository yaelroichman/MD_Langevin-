classdef additionalData < handle
    properties
        wallMoveSteps % The step numbers in which the wall was moved
        closestParticlePositions % The positions of the closest particle in each wall move
        closestParticleDistances % The distance ofthe closest particle from the right wall
        newWallPositions % The positions to which the wall was moved
        wallMoveInd;
        wallCheckInd;
    end
end