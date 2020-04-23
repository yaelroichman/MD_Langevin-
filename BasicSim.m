cfg = simConfig;
cfg.N = 1e8;
cfg.Dt = 1e-4;
cfg.sampleRate = 1;
cfg.numOfParticles = 2;
cfg.initPositions = [0 0; 1 1];
cfg.eta = 1;
cfg.R = 1;
cfg.T = 300;

cfg.useWalls = false;
cfg.useParticleRepulsion = false;
cfg.useHydro = false;
cfg.useTraps = false;

cfg.displayLive = true;

cfg.saveFoldername = 'a';
cfg.savePeriod = 1e5;

particlePositions = MDSim(cfg, @computeForces, @checkIfMoveWall, @moveWall, @printCurrStep, 1)


