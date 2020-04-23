function f = checkIfMoveWall(cfg, currStepData, addedData)
    rightmostPosition = max(currStepData.particlePositions(:,1)) + cfg.R(1);
    addedData.closestParticlePositions(addedData.wallCheckInd) = rightmostPosition;
    addedData.closestParticleDistances(addedData.wallCheckInd) = currStepData.wallPositionsX(2) - rightmostPosition;
    addedData.wallCheckInd = addedData.wallCheckInd + 1;
    f = false;
