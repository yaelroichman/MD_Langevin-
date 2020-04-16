function particlePositions = infoChamber(N,Dt,sampleRate,R,T,eta, lx, ly, numOfParticles, wallShrink,saveFoldername, displayLive)
    %% Setting the configuration for the simulation
    cfg = simConfig;
    cfg.N = N;
    cfg.Dt = Dt;
    cfg.sampleRate = sampleRate;
    cfg.numOfParticles = numOfParticles;
    cfg.eta = eta;
    cfg.R = ones(1,numOfParticles).*R;
    cfg.T = 300;
    cfg.wallPositionsX = [-lx/2, lx/2 - wallShrink];
    cfg.wallPositionsY = [-ly/2, ly/2];
    cfg.useWalls = true;
    cfg.useParticleRepulsion = true;
    cfg.WCAEpsilon = 0.2*1.38e-23*cfg.T; % 0.2 kBT seems to work well
    cfg.useHydro = true;
    cfg.displayLive = displayLive;
    cfg.saveFoldername = saveFoldername;
    cfg.savePeriod = 1e5;

    %% Creating the additional info object
    wallData = additionalData;
    maxWallMoves = 100;
    wallData.wallMoveSteps = zeros(1,maxWallMoves);
    wallData.newWallPositions = zeros(1,maxWallMoves);
    wallData.newWallPositions(1) = cfg.wallPositionsX(2);
    wallData.wallMoveSteps(1) = 0;
    wallData.wallMoveInd = 2;
    %% Randomizing particle starting positions
    rng(0,'twister');
    [particlesX, particlesY] = randomizePositions(cfg.wallPositionsX, cfg.wallPositionsY, numOfParticles, R);
    initPositions = zeros(2,numOfParticles);
    initPositions(1,:) = particlesX;
    initPositions(2,:) = particlesY;
    cfg.initPositions = initPositions;
    %% Running the simulation
    particlePositions = MDSim(cfg, @computeForces, @(x,y,z) false, @moveWall, @printCurrStep, wallData);
    %% Checking the mean densities
    pixelsPerLength = 10/1e-6
    densityHeatMap(particlePositions, cfg.R,...
                   cfg.wallPositionsX, cfg.wallPositionsY,...
                   pixelsPerLength, true, true, cfg.saveFoldername);
    %% Showing the results
    wallData.wallMoveSteps = wallData.wallMoveSteps(wallData.wallMoveSteps ~= 0);
    wallData.newWallPositions = wallData.newWallPositions(wallData.newWallPositions ~= 0);
    sampledWallMoves = wallData.wallMoveSteps / (sampleRate / Dt);
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
        plot([wallData.newWallPositions(movedInd),wallData.newWallPositions(movedInd)],[cfg.wallPositionsY(1),cfg.wallPositionsY(2)],'-');
        plot([cfg.wallPositionsX(1),wallData.newWallPositions(movedInd)],[cfg.wallPositionsY(1),cfg.wallPositionsY(1)],'-');
        plot([cfg.wallPositionsX(1),wallData.newWallPositions(movedInd)],[cfg.wallPositionsY(2),cfg.wallPositionsY(2)],'-');
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
%     % create the video writer with 1 fps
%     writerObj = VideoWriter(strcat(cfg.saveFoldername, '/movie.avi'));
%     writerObj.FrameRate = sampleRate;
%     % set the seconds per image
%     % open the video writer
%     open(writerObj);
%     % write the frames to the video
%     for i=1:length(F)
%         % convert the image to a frame
%         frame = F(i) ;    
%         writeVideo(writerObj, frame);
%     end
%     % close the writer object
%     close(writerObj);