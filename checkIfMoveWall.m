function f = checkIfMoveWall(cfg, currStepData, addedData)
    f = mod(currStepData.stepNum, 5e4) == 0;