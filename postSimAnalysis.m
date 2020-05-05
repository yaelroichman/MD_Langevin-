function postSimAnalysis(cfg, addedData)
    %% Reading the output file
    particlePositionsX = dlmread(strcat(cfg.saveFoldername, '/pos_x.csv'));
    particlePositionsY = dlmread(strcat(cfg.saveFoldername, '/pos_y.csv'));
    particlePositions = zeros(size(particlePositionsX, 1), cfg.numOfParticles, 2);
    particlePositions(:,:,2) = particlePositionsY;
    particlePositions(:,:,1) = particlePositionsX;
    %% Checking the percentages of movement options
%     addedData.closestParticlePositions = addedData.closestParticlePositions(~isnan(addedData.closestParticlePositions)); 
%     binWidth = 2e-7;
%     histEdges = 0:binWidth:max(addedData.closestParticleDistances);
%     figure(2);
%     histogram(addedData.closestParticleDistances, 'BinEdges', histEdges, 'Normalization', 'probability');
%     xlabel('Distance [m]');
%     ylabel('percentage of time in distance');
%     saveas(gcf, strcat(cfg.saveFoldername,'/distances.png'));
    %% Checking the mean densities
    pixelsPerLength = 10/1e-6;
    densityHeatMap(particlePositions, cfg.R,...
                   cfg.xlimits, cfg.ylimits,...
                   pixelsPerLength, false, false, cfg.saveFoldername);
    %% Showing the tracks
    addedData.wallMoveSteps = addedData.wallMoveSteps(~isnan(addedData.wallMoveSteps));
    addedData.newWallPositions = addedData.newWallPositions(~isnan(addedData.newWallPositions));
    sampledWallMoves = addedData.wallMoveSteps / (cfg.sampleRate / cfg.Dt);
    movedInd = 1;
    ColorSet = varycolor(cfg.numOfParticles);
    t = [0:cfg.Dt:(cfg.N-1)*cfg.Dt];
    figure(5)
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
    saveas(gcf, strcat(cfg.saveFoldername,'/tracks.png'));
    %% Saving the diffusion as a movie
%     numOfFrames = size(particlePositions,1);
%     for i = 1:numOfFrames
%         figure(1)
%         title(i);
%         plot([cfg.wallPositionsX(1),cfg.wallPositionsX(1)],[cfg.wallPositionsY(1),cfg.wallPositionsY(2)],'-');
%         hold on
%         plot([addedData.newWallPositions(movedInd),addedData.newWallPositions(movedInd)],[cfg.wallPositionsY(1),cfg.wallPositionsY(2)],'-');
%         plot([cfg.wallPositionsX(1),addedData.newWallPositions(movedInd)],[cfg.wallPositionsY(1),cfg.wallPositionsY(1)],'-');
%         plot([cfg.wallPositionsX(1),addedData.newWallPositions(movedInd)],[cfg.wallPositionsY(2),cfg.wallPositionsY(2)],'-');
%         viscircles([particlePositions(i,:,1)',...
%                     particlePositions(i,:,2)'],...
%             ones(cfg.numOfParticles,1).*cfg.R(1));
%         hold off
%           F(i) = getframe(gcf);
%           drawnow
%         if movedInd < length(sampledWallMoves) && mod(i,sampledWallMoves(movedInd)) == 0
%             movedInd = movedInd + 1;
%         end
%     end
%     saveMovie(F, 1/cfg.sampleRate, strcat(cfg.saveFoldername, '/movie.avi'))