function moveWall(cfg, currStepData, addedData)
    currWidth = abs(currStepData.wallPositionsX(2) - currStepData.wallPositionsX(1));
    newWidth = 0.9*currWidth;
    
    currStepData.wallPositionsX(2) = rightmostPosition;
    addedData.wallMoveSteps(addedData.wallMoveInd) = currStepData.stepNum;
    addedData.newWallPositions(addedData.wallMoveInd) = rightmostPosition;
    addedData.wallMoveInd = addedData.wallMoveInd + 1;