function printCurrStep(cfg, currStepData, addedData)
    if cfg.useWalls
        plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(1)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(2)],'-');
        hold on
        plot([currStepData.wallPositionsX(2),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(2)],'-');
        plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(1)],'-');
        plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(2),currStepData.wallPositionsY(2)],'-');
    end

    viscircles([currStepData.particlePositions(:,1),...
                currStepData.particlePositions(:,2)],...
        ones(cfg.numOfParticles,1).*cfg.R(1));
    hold off
    pause(0.001); 