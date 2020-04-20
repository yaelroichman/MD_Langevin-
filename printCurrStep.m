function printCurrStep(cfg, currStepData, addedData)
    cla     
    if cfg.useWalls
        plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(1)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(2)],'-');
        hold on
        plot([currStepData.wallPositionsX(2),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(2)],'-');
        plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(1)],'-');
        plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(2),currStepData.wallPositionsY(2)],'-');
    end
   xlim([-cfg.xlimit,cfg.xlimit])
   ylim([-cfg.ylimit,cfg.ylimit])
    viscircles([currStepData.particlePositions(:,1),...
                currStepData.particlePositions(:,2)],...
        ones(cfg.numOfParticles,1).*cfg.R(1));
  %  hold off
    pause(0.001); 