function particlePositions = MDSim(cfg, forcesFunc, feedbackCheckFunc, feedbackFunc, addedData)
close all
kB = 1.38e-23; % Boltzmann constant [J/K]
R = cfg.R(1); % Currently works only with one constant R
gamma = 6*pi*R*cfg.eta; % friction coefficient
D = kB*cfg.T/gamma; % diffusion coefficient
d = 2; % The dimension of the problem. Currently ONLY WORKS FOR 2d!!
epsilon = cfg.WCAEpsilon; % Lennard-Jones interaction strength
wallPositionsX = cfg.wallPositionsX;
wallPositionsY = cfg.wallPositionsY;
sampleRate = cfg.sampleRate; %fps
samplePeriod = round(sampleRate / cfg.Dt);
particlePositions = zeros(2,cfg.numOfParticles,cfg.N/samplePeriod);
particlePositions(:,:,1) = cfg.initPositions;
%% Plotting the initial placements
hold on
plot([wallPositionsX(1),wallPositionsX(1)],[wallPositionsY(1),wallPositionsY(2)],'-');
plot([wallPositionsX(2),wallPositionsX(2)],[wallPositionsY(1),wallPositionsY(2)],'-');
plot([wallPositionsX(1),wallPositionsX(2)],[wallPositionsY(1),wallPositionsY(1)],'-');
plot([wallPositionsX(1),wallPositionsX(2)],[wallPositionsY(2),wallPositionsY(2)],'-');
% plot(particlesX,particlesY,'o');
particlesX = squeeze(particlePositions(1,:,1))';
particlesY = squeeze(particlePositions(2,:,1))';
viscircles([particlesX, particlesY],...
            ones(cfg.numOfParticles,1).*R);
hold off
xlabel('x [m]');
ylabel('y [m]');
title('Initial placement');
pause(0.01);
%% Checking whether to run hydrosynamic interactions
if ~cfg.useHydro
    Dx = ones(1, numOfParticles).*D;
    Dy = Dx;
    Ax = ones(1, numOfParticles).*sqrt(2*D);
    Ay = Ax;
end
%% Running the simulation
sampleInd = 2;
printPeriod = 1e5;
currStepData = simStepData;
currStepData.particlePositions = (squeeze(particlePositions(:,:,1)));
currStepData.wallPositionsX = wallPositionsX;
currStepData.wallPositionsY = wallPositionsY;
for i = 2:1:cfg.N
    currStepData.stepNum = i;
    %% Printing data
    if mod(i, printPeriod) == 0
        i
    end
    if mod(i, samplePeriod) == 0
        if cfg.displayLive
            if cfg.useWalls
                plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(1)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(2)],'-');
                hold on
                plot([currStepData.wallPositionsX(2),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(2)],'-');
                plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(1)],'-');
                plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(2),currStepData.wallPositionsY(2)],'-');
            end

            viscircles([currStepData.particlePositions(1,:)',...
                        currStepData.particlePositions(2,:)'],...
                ones(cfg.numOfParticles,1).*R);
            hold off
            pause(0.001);            
        end

        particlePositions(:,:,sampleInd) = currStepData.particlePositions;
        sampleInd = sampleInd + 1;
    end
    %% Forces computation
    [fx, fy] = ...
        forcesFunc(cfg, currStepData, addedData);
    %% Hydrodynamic interactions computation
    if cfg.useHydro
        DMat = rotnePrager(currStepData.particlePositions',R, D);
        rootMat = chol(DMat,'lower');
        DVec = DMat*ones(cfg.numOfParticles.*d,1);
        AVec = rootMat*ones(cfg.numOfParticles.*d,1);
        Dx = DVec(1:2:end);
        Dy = DVec(2:2:end);
        Ax = AVec(1:2:end);
        Ay = AVec(2:2:end);
    end
    %% Running the step
    currStepData.particlePositions(1,:) = currStepData.particlePositions(1,:) +...
                    Ax'.*sqrt(cfg.Dt).*randn(1,cfg.numOfParticles) +...
                    (Dx'./(kB.*cfg.T)).*fx'.*cfg.Dt;
    currStepData.particlePositions(2,:) = currStepData.particlePositions(2,:) +...
                    Ay'.*sqrt(cfg.Dt).*randn(1,cfg.numOfParticles) +...
                    (Dy'./(kB.*cfg.T)).*fy'.*cfg.Dt;
    %% Checking whether to apply a feedback to the system, and applying it if necessary
    if feedbackCheckFunc(cfg, currStepData, addedData)
        feedbackFunc(cfg, currStepData, addedData);
    end
end