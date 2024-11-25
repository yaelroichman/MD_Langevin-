function f = checkIfMoveTraps(cfg, currStepData, addedData)
if mod(currStepData.stepNum,cfg.switchingPeriod)==0 
    f=true;
else 
    f=false;
end