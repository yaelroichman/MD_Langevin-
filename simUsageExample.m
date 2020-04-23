    %% Setting the configuration for the simulation
    cfg = simConfig; % This object is used to set the initial conditions of the simulation
    cfg.N = 1e5;
    cfg.Dt = 1e-4;
    cfg.sampleRate = 1/60;
    cfg.numOfParticles = 3;
    cfg.eta = 0.001;
    cfg.R = ones(1,cfg.numOfParticles).*R;
    cfg.T = 300;
    cfg.useWalls = false;
    cfg.useParticleRepulsion = true;
    cfg.WCAEpsilon = 0.2*1.38e-23*cfg.T; % 0.2 kBT seems to work well
    cfg.useHydro = true;
    cfg.displayLive = true;
    cfg.saveFoldername = 'MDSim Example 1';
    cfg.savePeriod = 1e3;
    cfg.initPositions = 

    %% Creating the additional info object
    addedData = additionalData;
    % This object is an output parameter used for tracking additional
    % elements of the simulation (other than the particle tracks)
    %% Randomizing particle starting positions
    rng('shuffle');
    [particlesX, particlesY] = randomizePositions(cfg.wallPositionsX, cfg.wallPositionsY, numOfParticles, R);
    initPositions = zeros(2,numOfParticles);
    initPositions(1,:) = particlesX;
    initPositions(2,:) = particlesY;
    cfg.initPositions = initPositions;
    %% Running the simulation
    particlePositions = MDSim(cfg, @computeForces, @(x,y,z) false, @moveWall, @printCurrStep, addedData);
    %% Checking the mean densities
    pixelsPerLength = 10/1e-6;
    densityHeatMap(particlePositions, cfg.R,...
                   cfg.wallPositionsX, cfg.wallPositionsY,...
                   pixelsPerLength, true, true, cfg.saveFoldername);
    %% Showing the results
    addedData.wallMoveSteps = addedData.wallMoveSteps(~isnan(addedData.wallMoveSteps));
    addedData.newWallPositions = addedData.newWallPositions(~isnan(addedData.newWallPositions));
    sampledWallMoves = addedData.wallMoveSteps / (sampleRate / Dt);
    movedInd = 1;
    ColorSet = varycolor(cfg.numOfParticles);
    t = [0:cfg.Dt:(N-1)*cfg.Dt];
    figure
    hold on
    plot([cfg.wallPositionsX(1),cfg.wallPositionsX(1)],[cfg.wallPositionsY(1),cfg.wallPositionsY(2)],'-');
    plot([cfg.wallPositionsX(2),cfg.wallPositionsX(2)],[cfg.wallPositionsY(1),cfg.wallPositionsY(2)],'-');
    plot([cfg.wallPositionsX(1),cfg.wallPositionsX(2)],[cfg.wallPositionsY(1),cfg.wallPositionsY(1)],'-');
    plot([cfg.wallPositionsX(1),cfg.wallPositionsX(2)],[cfg.wallPositionsY(2),cfg.wallPositionsY(2)],'-');
    for i=1:cfg.numOfParticles
        currParticleTrackX = squeeze(particlePositions(:,i,1));
        currParticleTrackY = squeeze(particlePositions(:,i,2));
        plot(currParticleTrackX, currParticleTrackY, '.', 'Color', ColorSet(i,:));
    end

    xlabel('x [m]')
    ylabel('y [m]')
    title('tracks')

    %% Saving the diffusion as a movie
    numOfFrames = size(particlePositions,1);
    for i = 1:numOfFrames

        figure(1)  
        plot([cfg.wallPositionsX(1),cfg.wallPositionsX(1)],[cfg.wallPositionsY(1),cfg.wallPositionsY(2)],'-');
        hold on
        plot([addedData.newWallPositions(movedInd),addedData.newWallPositions(movedInd)],[cfg.wallPositionsY(1),cfg.wallPositionsY(2)],'-');
        plot([cfg.wallPositionsX(1),addedData.newWallPositions(movedInd)],[cfg.wallPositionsY(1),cfg.wallPositionsY(1)],'-');
        plot([cfg.wallPositionsX(1),addedData.newWallPositions(movedInd)],[cfg.wallPositionsY(2),cfg.wallPositionsY(2)],'-');
        viscircles([particlePositions(i,:,1)',...
                    particlePositions(i,:,2)'],...
            ones(numOfParticles,1).*cfg.R(1));
        hold off
          F(i) = getframe(gcf);
          drawnow
        if movedInd < length(sampledWallMoves) && mod(i,sampledWallMoves(movedInd)) == 0
            movedInd = movedInd + 1;
        end
    end
    saveMovie(F, 1/cfg.sampleRate, strcat(cfg.saveFoldername, '/movie.avi'))