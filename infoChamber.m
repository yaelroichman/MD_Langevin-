function particlePositions = infoChamber(N,Dt,sampleRate,R,T,eta, lx, ly, numOfParticles, wallShrink,saveFoldername, displayLive)
    %% Setting the configuration for the simulation
    cfg = simConfig;
    cfg.N = N;
    cfg.Dt = Dt;
    cfg.sampleRate = sampleRate;
    cfg.numOfParticles = numOfParticles;
    cfg.eta = eta;
    cfg.R = ones(1,numOfParticles).*R;
    cfg.T = T;
    cfg.wallPositionsX = [-lx/2, lx/2 - wallShrink];
    cfg.wallPositionsY = [-ly/2, ly/2];
    cfg.xlimits=cfg.wallPositionsX.*1.05;
    cfg.ylimits=cfg.wallPositionsY*1.05;
    cfg.useWalls = true;
    cfg.useParticleRepulsion = true;
    cfg.WCAEpsilon = 0.2*physconst('boltzmann')*cfg.T; % 0.2 kBT seems to work well
    cfg.useHydro = true;
    cfg.displayLive = displayLive;
    cfg.saveFoldername = saveFoldername;
    cfg.savePeriod = cfg.N ./ 10; %save to disk

    %% Creating the additional info object
    wallData = additionalData;
    maxWallMoves = 100;
    wallData.wallMoveSteps = nan(1,maxWallMoves);
    wallData.newWallPositions = nan(1,maxWallMoves);
    wallData.newWallPositions(1) = cfg.wallPositionsX(2);
    wallData.wallMoveSteps(1) = 0;
    wallData.wallMoveInd = 2;
    wallData.wallCheckInd = 1;
    wallData.closestParticlePositions = nan(cfg.N,1);
    wallData.closestParticleDistances = nan(cfg.N,1);
    %% Randomizing particle starting positions
    rng('shuffle');
    [particlesX, particlesY] = randomizePositions(cfg.wallPositionsX, cfg.wallPositionsY, numOfParticles, R);
    initPositions = zeros(numOfParticles, 2);
    initPositions(:,1) = particlesX;
    initPositions(:,2) = particlesY;
    cfg.initPositions = initPositions;
    %% Running the simulation
    particlePositions = MDSim(cfg, @computeForces, @checkIfMoveWall, @moveWall, @printCurrStep, wallData);
    
    %% Running the post-simulation analysis
    postSimAnalysis(cfg, wallData);