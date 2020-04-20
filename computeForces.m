function [fx, fy] = computeForces(cfg, currStepData, addedData)

fx = 0;
fy = 0;

if cfg.useWalls || cfg.useParticleRepulsion
    [fxWCA, fyWCA] = getWCAForces(currStepData.particlePositions, cfg.R(1), cfg.WCAEpsilon,...
        cfg.useParticleRepulsion, cfg.useWalls,...
        currStepData.wallPositionsX, currStepData.wallPositionsY);
    fx = fx + fxWCA;
    fy = fy + fyWCA;
end

if cfg.useTraps
    [fxTraps, fyTraps] = getTrapForces(currStepData.particlePositions, currStepData.trapPositions, cfg.A, cfg.s);
    fx = fx + fxTraps;
    fy = fy + fyTraps;
end